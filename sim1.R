#This script does this this this....
#conda activate env1
(base) usman@bioinfo:~$ pwd
/home/usman
(base) usman@bioinfo:~$ mv Files/* /data/sim_breeding/data/
  (base) usman@bioinfo:~$ ls Files/
  (base) usman@bioinfo:~$ ls /data/sim_breeding/
  code  data
(base) usman@bioinfo:~$ ls /data/sim_breeding/data/
  accessions.rds          ath_all_new.fam               ath_all_new_maf_ldpruned.bim  MAC_df.rds                  Pheno.rds
all_nucl_genes_bed.rds  ath_all_new_maf.bed           ath_all_new_maf_ldpruned.fam  MAC_matrix_with_header.rds
ath_all_new.bim         ath_all_new_maf_ldpruned.bed  ath_all_new_maf_ldpruned.map  Pheno_pla.Rdata
(base) usman@bioinfo:~$ conda activate env1
(env1) usman@bioinfo:~$ R
######################################################################################################################################################
setwd('/data/sim_breeding/code')
data_path="/data/sim_breeding/data/"
######################################################################################################################################################
library(AlphaSimR)
library(data.table)
library(BEDMatrix)
######################################################################################################################################################
#load genotype data
paste0(data_path,'ath_all_new_maf_ldpruned')
data_path="/data/sim_breeding/data/"

geno=BEDMatrix(paste0(data_path,'ath_all_new_maf_ldpruned'),simple_names=TRUE)

dim(geno)

geno[1:5,1:5]


######################################################################################################################################################
> founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)

pop = newPop(founderPop, simParam=SP)
aa(pop, simParam=SP)
Trait1
[1,]      0
[2,]      0
[3,]      0
[4,]      0
[5,]      0
[6,]      0
[7,]      0
[8,]      0
[9,]      0
[10,]     0

######################################################################################################################################################
> founderPop = quickHaplo(nInd=10,nChr=1,segSites=2)

haplo = matrix(sample(x=0:1, size=20, replace=TRUE), ncol=1)
founderPop2 = addSegSite(founderPop, siteName="x", chr=1, mapPos=0.5, haplo=haplo)
pullSegSiteHaplo(founderPop2)
1_1 x 1_2
1_1    1 1   1
1_2    0 1   1
2_1    0 1   1
2_2    1 1   0
3_1    0 0   0
3_2    1 0   1
4_1    1 1   0
4_2    1 0   0
5_1    0 1   0
5_2    1 1   1
6_1    1 0   1
6_2    0 1   0
7_1    1 1   0
7_2    0 1   0
8_1    0 1   1
8_2    0 1   0
9_1    1 0   1
9_2    1 0   1
10_1   0 0   1
10_2   1 1   0

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)
#Create population
pop = newPop(founderPop, simParam=SP)
bv(pop, simParam=SP)
Trait1
[1,] -0.69381762
[2,] -0.17407552
[3,] -1.43636958
[4,]  1.63378298
[5,]  1.48207197
[6,]  0.31952742
[7,] -1.53934734
[8,]  0.26812380
[9,]  0.01537012
[10,]  0.12473378

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)
#Create population
pop = newPop(founderPop, simParam=SP)
dd(pop, simParam=SP)
Trait1
[1,]  0.37155653
[2,]  0.40823264
[3,] -0.27997318
[4,] -0.20817194
[5,] -0.29072601
[6,] -0.14532698
[7,]  0.14450819
[8,]  0.11812909
[9,] -0.17612079
[10,]  0.05789246

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)
#Create population
pop = newPop(founderPop, simParam=SP)
pop@ebv = matrix(rnorm(pop@nInd), nrow=pop@nInd, ncol=1)
ebv(pop)
[,1]
[1,] -0.09461147
[2,] -1.87477769
[3,]  0.59546820
[4,] -3.49007480
[5,] -0.12713337
[6,] -0.71055880
[7,] -0.01359387
[8,]  0.11726226
[9,]  1.05183897
[10,] -2.20098403

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)
SP$addSnpChip(10)
#Create population
pop = newPop(founderPop, simParam=SP)
#Run GS model and set EBV
ans = fastRRBLUP(pop, simParam=SP)
pop = setEBV(pop, ans, simParam=SP)
#Evaluate accuracy
cor(gv(pop), ebv(pop))
est_GV_Trait1
Trait1     0.5393535

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
getGenMap(founderPop)
id chr       pos
1   1_1   1 0.0000000
2   1_2   1 0.1111111
3   1_3   1 0.2222222
4   1_4   1 0.3333333
5   1_5   1 0.4444444
6   1_6   1 0.5555556
7   1_7   1 0.6666667
8   1_8   1 0.7777778
9   1_9   1 0.8888889
10 1_10   1 1.0000000

######################################################################################################################################################
> # Create a founder population
  founderPop = quickHaplo(2,1,2)
# Set simulation parameters
SP = SimParam$new(founderPop)
# Create a population
pop = newPop(founderPop, simParam=SP)
# Get the pedigree
getPed(pop)
# Returns NULL when a population lacks a pedigree
getPed(founderPop)
id mother father
1  1      0      0
2  2      0      0
NULL

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(5)
#Pull SNP map
getQtlMap(trait=1, simParam=SP)
id chr site       pos
1_1   1_1   1    1 0.0000000
1_6   1_6   1    6 0.5555556
1_7   1_7   1    7 0.6666667
1_9   1_9   1    9 0.8888889
1_10 1_10   1   10 1.0000000

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addSnpChip(5)
#Pull SNP map
getSnpMap(snpChip=1, simParam=SP)
id chr site       pos
1_2   1_2   1    2 0.1111111
1_3   1_3   1    3 0.2222222

1_4   1_4   1    4 0.3333333
1_6   1_6   1    6 0.5555556
1_10 1_10   1   10 1.0000000

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)
#Create population
pop = newPop(founderPop, simParam=SP)
gv(pop)
Trait1
[1,] -2.0465832
[2,] -0.2637463
[3,]  0.9438442
[4,]  0.8874569
[5,]  0.2305051
[6,]  0.4848081
[7,]  1.4199048
[8,] -2.0549617
[9,]  0.2835943
[10,]  0.1151776

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)
#Create population
pop = newPop(founderPop, simParam=SP)
meanG(pop)
Trait1
1.620926e-15

######################################################################################################################################################
> #Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)
#Create population
pop = newPop(founderPop, simParam=SP)
meanP(pop)
Trait1
0.03145706

######################################################################################################################################################
>#Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)
SP$addSnpChip(10)
#Create population
pop = newPop(founderPop, simParam=SP)
#Run GS model and set EBV
ans = RRBLUP(pop, simParam=SP)
pop = setEBV(pop, ans, simParam=SP)
#Evaluate accuracy
cor(gv(pop), ebv(pop))
est_GV_Trait1
Trait1     0.7661872

######################################################################################################################################################
>#Create founder haplotypes
  founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.8)
#Create population
pop = newPop(founderPop, simParam=SP)

#Select top 5 (directional selection)
pop2 = selectInd(pop, 5, simParam=SP)
hist(pop@pheno); abline(v = pop@pheno, lwd = 2)
abline(v = pop2@pheno, col = "red", lwd = 2)

#Select 5 most deviating from an optima (disruptive selection)
squaredDeviation = function(x, optima = 0) (x - optima)^2
pop3 = selectInd(pop, 5, simParam=SP, trait = squaredDeviation, selectTop = TRUE)
hist(pop@pheno); abline(v = pop@pheno, lwd = 2)
abline(v = pop3@pheno, col = "red", lwd = 2)

#Select 5 least deviating from an optima (stabilising selection)
pop4 = selectInd(pop, 5, simParam=SP, trait = squaredDeviation, selectTop = FALSE)
hist(pop@pheno); abline(v = pop@pheno, lwd = 2)
abline(v = pop4@pheno, col = "red", lwd = 2)

#####################################################################################################################################################
#Create founder Population
founderPop = qickHaplo(nInd=50, nChr=5, segSites=100, inbred=TRUE)

#Select Parameters
SP = SimParam$new(founderPop)
SP$addTraitA(100)
SP$setVarE(h2=0.5)

#Select Parental lines
Parents=newPop(founderPop)
F1=randCross(Parents, 200)

#In the second and third years, the DH lines are produced
HDRW=makeDH(F1, 100)
HDRW=setPheno(HDRW, H2=0.1)

#In the fourth year, the best entries in the HDRW nursery are selected and evaluated in a preliminary yield trial (PYT).
PYT=selectWithinFam(HDRW, 10)
PYT=setPheno(PYT)

#In the fifth year, the best PYT entries are selected and evaluated in an advanced yield trial (AYT).
AYT=selectInd(PYT, 100)
AYT=setPheno(AYT, reps = 4)

#In the sixth year, the best AYT entries are selected and evaluated in an elite yield trial (EYT).
EYT=selectInd(AYT, 20)
EYT=setPheno(EYT, reps = 16)

#In the seventh year, the best-performing EYT entry is chosen for release as a variety.
Variety=selectInd(EYT, 2)

#The final step is to evaluate the simulation results.
yield=list(Parents=gv(Parents), F1=gv(F1),HDRW=gv(HDRW),PYT=gv(PYT),AYT=gv(AYT),EYT=gv(EYT),Variety=gv(Variety))
boxplot(yield, ylab="Genetic value")

######################################################################################################################################################
#Create founder Population
founderPop = runMacs(nInd=100, nChr=5, segSites=1000, inbred=TRUE)

#Select Parameters
SP = SimParam$new(founderPop)
SP$addTraitA(1000)
SP$setVarE(h2=0.5)

#Select Parental lines
Parents=newPop(founderPop)
F1 = selectCross(Parents, nInd=10, nCrosses=100, simParam=SP)

#In the second and third years, the DH lines are produced
HDRW=makeDH(F1, 50)
HDRW=setPheno(HDRW, H2=0.1)

#In the fourth year, the best entries in the HDRW nursery are selected and evaluated in a preliminary yield trial (PYT).
PYT=selectWithinFam(HDRW, 10)
PYT=setPheno(PYT)

#In the fifth year, the best PYT entries are selected and evaluated in an advanced yield trial (AYT).
AYT=selectInd(PYT, 50)
AYT=setPheno(AYT, reps = 5)

#In the sixth year, the best AYT entries are selected and evaluated in an elite yield trial (EYT).
EYT=selectInd(AYT, 25)
EYT=setPheno(EYT, reps = 10)

#In the seventh year, the best-performing EYT entry is chosen for release as a variety.
Variety=selectInd(EYT, 5)

#The final step is to evaluate the simulation results.
yield=list(Parents=gv(Parents), F1=gv(F1),HDRW=gv(HDRW),PYT=gv(PYT),AYT=gv(AYT),EYT=gv(EYT),Variety=gv(Variety))
boxplot(yield, ylab="Genetic value")

######################################################################################################################################################
#Create founder Population
founderPop = runMacs(nInd=100, nChr=5, segSites=1000, inbred=TRUE)

#Select Parameters
SP = SimParam$new(founderPop)
SP$addTraitA(1000)
SP$setVarE(h2=0.2)

#Select Parental lines
Parents=newPop(founderPop)
F1 = hybridCross(Parents,Parents, simParam=SP)

#In the second and third years, the DH lines are produced
HDRW=makeDH(F1, 50)
HDRW=setPheno(HDRW, H2=0.1)

#In the fourth year, the best entries in the HDRW nursery are selected and evaluated in a preliminary yield trial (PYT).
PYT=selectWithinFam(HDRW, 10)
PYT=setPheno(PYT)

#In the fifth year, the best PYT entries are selected and evaluated in an advanced yield trial (AYT).
AYT=selectInd(PYT, 50)
AYT=setPheno(AYT, reps = 5)

#In the sixth year, the best AYT entries are selected and evaluated in an elite yield trial (EYT).
EYT=selectInd(AYT, 25)
EYT=setPheno(EYT, reps = 10)

#In the seventh year, the best-performing EYT entry is chosen for release as a variety.
EYT=setPheno(EYT, reps = 10)

#The final step is to evaluate the simulation results.
variety=selectInd(EYT,1)
yield=list(Parents=gv(Parents), F1=gv(F1),HDRW=gv(HDRW),PYT=gv(PYT),AYT=gv(AYT),EYT=gv(EYT),Variety=gv(Variety))
boxplot(yield, ylab="Genetic value")

######################################################################################################################################################
#Create founder Population
founderPop = runMacs(nInd=10, nChr=5, segSites=100, inbred=TRUE)

#Select Parameters
SP = SimParam$new(founderPop)
SP$addTraitA(100)
SP$setVarE(h2=0.5)

#Select Parental lines
Parents=newPop(founderPop)
F1=randCross(Parents, 200)

# Assuming 'Parents' is your simulated population
Parents=setPheno(Parents)

#Get QTL map
getQtlMap(trait=1, simParam=SP)
id chr site          pos
1_1     1_1   1    1 0.0000000000
1_2     1_2   1    2 0.0020858806
1_3     1_3   1    3 0.0127535685
1_4     1_4   1    4 0.0178394750
1_5     1_5   1    5 0.0185640522
1_6     1_6   1    6 0.0268082648
1_7     1_7   1    7 0.0605597389
1_8     1_8   1    8 0.0693386949
1_9     1_9   1    9 0.0757111493
1_10   1_10   1   10 0.0808383930
1_11   1_11   1   11 0.0828605814
1_12   1_12   1   12 0.0865037303
1_13   1_13   1   13 0.0872835257
1_14   1_14   1   14 0.0982238514
1_15   1_15   1   15 0.0990340593
1_16   1_16   1   16 0.1028902172
1_17   1_17   1   17 0.1048074774
1_18   1_18   1   18 0.1129770262
1_19   1_19   1   19 0.1138824375
1_20   1_20   1   20 0.1143824382
1_21   1_21   1   21 0.1153014751
1_22   1_22   1   22 0.1168543223
1_23   1_23   1   23 0.1339037462
1_24   1_24   1   24 0.1342980575
1_25   1_25   1   25 0.1420993285
1_26   1_26   1   26 0.1706222647
1_27   1_27   1   27 0.1756959300
1_28   1_28   1   28 0.1769573106
1_29   1_29   1   29 0.1861900020
1_30   1_30   1   30 0.1890496711
1_31   1_31   1   31 0.1915184447
1_32   1_32   1   32 0.1916636404
1_33   1_33   1   33 0.1944367709
1_34   1_34   1   34 0.1992578342
1_35   1_35   1   35 0.2085820516
1_36   1_36   1   36 0.2209970370
1_37   1_37   1   37 0.2218447049
1_38   1_38   1   38 0.2233413193
1_39   1_39   1   39 0.2267242959
1_40   1_40   1   40 0.2369220953
1_41   1_41   1   41 0.2375942317
1_42   1_42   1   42 0.2494021901
1_43   1_43   1   43 0.2549979027
1_44   1_44   1   44 0.2841563231
1_45   1_45   1   45 0.2895993842
1_46   1_46   1   46 0.2913824396
1_47   1_47   1   47 0.3186561189
1_48   1_48   1   48 0.3265213033
1_49   1_49   1   49 0.3315902388
1_50   1_50   1   50 0.3326144422
1_51   1_51   1   51 0.3505760927
1_52   1_52   1   52 0.3538986202
1_53   1_53   1   53 0.3997174870
1_54   1_54   1   54 0.3999150394
1_55   1_55   1   55 0.3999324807
1_56   1_56   1   56 0.4037990219
1_57   1_57   1   57 0.4106504409
1_58   1_58   1   58 0.4110517807
1_59   1_59   1   59 0.4148658443
1_60   1_60   1   60 0.4465916731
1_61   1_61   1   61 0.4576013663
1_62   1_62   1   62 0.4678174073
1_63   1_63   1   63 0.4733643174
1_64   1_64   1   64 0.4938391324
1_65   1_65   1   65 0.5020613750
1_66   1_66   1   66 0.5275654034
1_67   1_67   1   67 0.5468278268
1_68   1_68   1   68 0.5495126471
1_69   1_69   1   69 0.5549575610
1_70   1_70   1   70 0.5549867828
1_71   1_71   1   71 0.6195886098
1_72   1_72   1   72 0.6651661232
1_73   1_73   1   73 0.6799439181
1_74   1_74   1   74 0.6955096952
1_75   1_75   1   75 0.6998540748
1_76   1_76   1   76 0.7058261550
1_77   1_77   1   77 0.7196051696
1_78   1_78   1   78 0.7370028699
1_79   1_79   1   79 0.7392703112
1_80   1_80   1   80 0.7561875878
1_81   1_81   1   81 0.7647596986
1_82   1_82   1   82 0.7719808536
1_83   1_83   1   83 0.7932331825
1_84   1_84   1   84 0.8039763492
1_85   1_85   1   85 0.8114305038
1_86   1_86   1   86 0.8407696220
1_87   1_87   1   87 0.8589327872
1_88   1_88   1   88 0.8907307609
1_89   1_89   1   89 0.8913760642
1_90   1_90   1   90 0.9106016919
1_91   1_91   1   91 0.9170276120
1_92   1_92   1   92 0.9178559413
1_93   1_93   1   93 0.9448661742
1_94   1_94   1   94 0.9454147566
1_95   1_95   1   95 0.9460343660
1_96   1_96   1   96 0.9579217999
1_97   1_97   1   97 0.9590253450
1_98   1_98   1   98 0.9616449795
1_99   1_99   1   99 0.9618689068
1_100 1_100   1  100 0.9652253942
2_1     2_1   2    1 0.0000000000
2_2     2_2   2    2 0.0006182975
2_3     2_3   2    3 0.0034498243
2_4     2_4   2    4 0.0243172489
2_5     2_5   2    5 0.0246390676
2_6     2_6   2    6 0.0834594439
2_7     2_7   2    7 0.0847419840
2_8     2_8   2    8 0.0854756914
2_9     2_9   2    9 0.0870764197
2_10   2_10   2   10 0.1064336158
2_11   2_11   2   11 0.1130979266
2_12   2_12   2   12 0.1132903270
2_13   2_13   2   13 0.1169571146
2_14   2_14   2   14 0.1457453774
2_15   2_15   2   15 0.1466000984
2_16   2_16   2   16 0.1516086528
2_17   2_17   2   17 0.1755659481
2_18   2_18   2   18 0.1778058947
2_19   2_19   2   19 0.1841740322
2_20   2_20   2   20 0.1855154199
2_21   2_21   2   21 0.2000268973
2_22   2_22   2   22 0.2019592613
2_23   2_23   2   23 0.2367504194
2_24   2_24   2   24 0.2374572970
2_25   2_25   2   25 0.2413602977
2_26   2_26   2   26 0.2424027189
2_27   2_27   2   27 0.2545067166
2_28   2_28   2   28 0.2656210125
2_29   2_29   2   29 0.2889461904
2_30   2_30   2   30 0.2899260683
2_31   2_31   2   31 0.2912400966
2_32   2_32   2   32 0.3018312606
2_33   2_33   2   33 0.3595964527
2_34   2_34   2   34 0.3657095400
2_35   2_35   2   35 0.3695243805
2_36   2_36   2   36 0.3751350405
2_37   2_37   2   37 0.3838417184
2_38   2_38   2   38 0.3845485702
2_39   2_39   2   39 0.3884030175
2_40   2_40   2   40 0.3900072835
2_41   2_41   2   41 0.3975598055
2_42   2_42   2   42 0.4006004929
2_43   2_43   2   43 0.4023503500
2_44   2_44   2   44 0.4038240663
2_45   2_45   2   45 0.4056269160
2_46   2_46   2   46 0.4081207513
2_47   2_47   2   47 0.4128268267
2_48   2_48   2   48 0.4225682287
2_49   2_49   2   49 0.4248997094
2_50   2_50   2   50 0.4267566587
2_51   2_51   2   51 0.4334730491
2_52   2_52   2   52 0.4705987635
2_53   2_53   2   53 0.5165897495
2_54   2_54   2   54 0.5226741358
2_55   2_55   2   55 0.5320306129
2_56   2_56   2   56 0.5522349634
2_57   2_57   2   57 0.5524332178
2_58   2_58   2   58 0.5578804773
2_59   2_59   2   59 0.5598790226
2_60   2_60   2   60 0.5628111546
2_61   2_61   2   61 0.5656520101
2_62   2_62   2   62 0.5878872223
2_63   2_63   2   63 0.5880761799
2_64   2_64   2   64 0.5888264120
2_65   2_65   2   65 0.6019877263
2_66   2_66   2   66 0.6032944528
2_67   2_67   2   67 0.6053858046
2_68   2_68   2   68 0.6161977683
2_69   2_69   2   69 0.6237539999
2_70   2_70   2   70 0.6273650850
2_71   2_71   2   71 0.6281438799
2_72   2_72   2   72 0.6667654774
2_73   2_73   2   73 0.6769672180
2_74   2_74   2   74 0.6830050971
2_75   2_75   2   75 0.6932049346
2_76   2_76   2   76 0.6967310667
2_77   2_77   2   77 0.6985296915
2_78   2_78   2   78 0.7179074828
2_79   2_79   2   79 0.7361574172
2_80   2_80   2   80 0.7500328636
2_81   2_81   2   81 0.7518450971
2_82   2_82   2   82 0.7634092886
2_83   2_83   2   83 0.7634547501
2_84   2_84   2   84 0.7652696571
2_85   2_85   2   85 0.7674729662
2_86   2_86   2   86 0.7730579977
2_87   2_87   2   87 0.7883580737
2_88   2_88   2   88 0.7977063499
2_89   2_89   2   89 0.8164550022
2_90   2_90   2   90 0.8298284336
2_91   2_91   2   91 0.8359369828
2_92   2_92   2   92 0.8408575787
2_93   2_93   2   93 0.8410144694
2_94   2_94   2   94 0.8558973505
2_95   2_95   2   95 0.9095538077
2_96   2_96   2   96 0.9384886566
2_97   2_97   2   97 0.9394027924
2_98   2_98   2   98 0.9658095742
2_99   2_99   2   99 0.9805932667
2_100 2_100   2  100 0.9872959429
3_1     3_1   3    1 0.0000000000
3_2     3_2   3    2 0.0077449280
3_3     3_3   3    3 0.0190329816
3_4     3_4   3    4 0.0216136905
3_5     3_5   3    5 0.0262208769
3_6     3_6   3    6 0.0265145620
3_7     3_7   3    7 0.0284723378
3_8     3_8   3    8 0.0288107389
3_9     3_9   3    9 0.0586002490
3_10   3_10   3   10 0.0946656947
3_11   3_11   3   11 0.1039884894
3_12   3_12   3   12 0.1116285761
3_13   3_13   3   13 0.1193483393
3_14   3_14   3   14 0.1291852187
3_15   3_15   3   15 0.1389310120
3_16   3_16   3   16 0.1400665340
3_17   3_17   3   17 0.1483691144
3_18   3_18   3   18 0.1515699409
3_19   3_19   3   19 0.1808878021
3_20   3_20   3   20 0.2125844182
3_21   3_21   3   21 0.2264005066
3_22   3_22   3   22 0.2279906762
3_23   3_23   3   23 0.2338852255
3_24   3_24   3   24 0.2512469262
3_25   3_25   3   25 0.2587449564
3_26   3_26   3   26 0.2599464973
3_27   3_27   3   27 0.2692450826
3_28   3_28   3   28 0.2722433677
3_29   3_29   3   29 0.2750455251
3_30   3_30   3   30 0.2866812056
3_31   3_31   3   31 0.2869211840
3_32   3_32   3   32 0.2907661394
3_33   3_33   3   33 0.2946510002
3_34   3_34   3   34 0.3100936452
3_35   3_35   3   35 0.3116332942
3_36   3_36   3   36 0.3308016992
3_37   3_37   3   37 0.3717278327
3_38   3_38   3   38 0.3765585126
3_39   3_39   3   39 0.3879926911
3_40   3_40   3   40 0.3940825588
3_41   3_41   3   41 0.3966807710
3_42   3_42   3   42 0.4311182302
3_43   3_43   3   43 0.4776335576
3_44   3_44   3   44 0.4804452401
3_45   3_45   3   45 0.4848220443
3_46   3_46   3   46 0.4866291277
3_47   3_47   3   47 0.4866935391
3_48   3_48   3   48 0.4869737567
3_49   3_49   3   49 0.4875819041
3_50   3_50   3   50 0.4913789178
3_51   3_51   3   51 0.5191910594
3_52   3_52   3   52 0.5222544812
3_53   3_53   3   53 0.5396732159
3_54   3_54   3   54 0.5452249825
3_55   3_55   3   55 0.5808912609
3_56   3_56   3   56 0.5957017804
3_57   3_57   3   57 0.6025531458
3_58   3_58   3   58 0.6243635420
3_59   3_59   3   59 0.6295527381
3_60   3_60   3   60 0.6297546239
3_61   3_61   3   61 0.6518140568
3_62   3_62   3   62 0.6618195925
3_63   3_63   3   63 0.6719437299
3_64   3_64   3   64 0.6730687116
3_65   3_65   3   65 0.6786247376
3_66   3_66   3   66 0.6831770356
3_67   3_67   3   67 0.6909991437
3_68   3_68   3   68 0.7240428649
3_69   3_69   3   69 0.7404326313
3_70   3_70   3   70 0.7429076368
3_71   3_71   3   71 0.7443482732
3_72   3_72   3   72 0.7559292195
3_73   3_73   3   73 0.7590916397
3_74   3_74   3   74 0.7648347647
3_75   3_75   3   75 0.7782912959
3_76   3_76   3   76 0.7846570197
3_77   3_77   3   77 0.7912099879
3_78   3_78   3   78 0.8393889815
3_79   3_79   3   79 0.8453947979
3_80   3_80   3   80 0.8456730610
3_81   3_81   3   81 0.8484916549
3_82   3_82   3   82 0.8522518947
3_83   3_83   3   83 0.8688758665
3_84   3_84   3   84 0.8764865792
3_85   3_85   3   85 0.8899376494
3_86   3_86   3   86 0.8924897761
3_87   3_87   3   87 0.9250017295
3_88   3_88   3   88 0.9256997000
3_89   3_89   3   89 0.9295578349
3_90   3_90   3   90 0.9296598421
3_91   3_91   3   91 0.9521134038
3_92   3_92   3   92 0.9550471619
3_93   3_93   3   93 0.9556593328
3_94   3_94   3   94 0.9570469358
3_95   3_95   3   95 0.9679786227
3_96   3_96   3   96 0.9759026904
3_97   3_97   3   97 0.9782567610
3_98   3_98   3   98 0.9793822402
3_99   3_99   3   99 0.9804532930
3_100 3_100   3  100 0.9835769500
4_1     4_1   4    1 0.0000000000
4_2     4_2   4    2 0.0151932530
4_3     4_3   4    3 0.0219118982
4_4     4_4   4    4 0.0235798848
4_5     4_5   4    5 0.0341098397
4_6     4_6   4    6 0.0341230189
4_7     4_7   4    7 0.0367879561
4_8     4_8   4    8 0.0404762620
4_9     4_9   4    9 0.0423487791
4_10   4_10   4   10 0.0465145227
4_11   4_11   4   11 0.0526193181
4_12   4_12   4   12 0.0533462736
4_13   4_13   4   13 0.0807490122
4_14   4_14   4   14 0.1006662857
4_15   4_15   4   15 0.1039203052
4_16   4_16   4   16 0.1071959180
4_17   4_17   4   17 0.1137069405
4_18   4_18   4   18 0.1341542400
4_19   4_19   4   19 0.1547881242
4_20   4_20   4   20 0.2384394580
4_21   4_21   4   21 0.2533559273
4_22   4_22   4   22 0.2534779159
4_23   4_23   4   23 0.2540957754
4_24   4_24   4   24 0.2659875829
4_25   4_25   4   25 0.2730311284
4_26   4_26   4   26 0.2933515090
4_27   4_27   4   27 0.2997083659
4_28   4_28   4   28 0.3063273723
4_29   4_29   4   29 0.3072762710
4_30   4_30   4   30 0.3360959310
4_31   4_31   4   31 0.3373017052
4_32   4_32   4   32 0.3384729067
4_33   4_33   4   33 0.3486042127
4_34   4_34   4   34 0.3547718616
4_35   4_35   4   35 0.3668972784
4_36   4_36   4   36 0.3816859161
4_37   4_37   4   37 0.3818688361
4_38   4_38   4   38 0.3871898957
4_39   4_39   4   39 0.3906512528
4_40   4_40   4   40 0.3925666720
4_41   4_41   4   41 0.4023596908
4_42   4_42   4   42 0.4235072800
4_43   4_43   4   43 0.4426600499
4_44   4_44   4   44 0.4656662633
4_45   4_45   4   45 0.4939207958
4_46   4_46   4   46 0.5063552133
4_47   4_47   4   47 0.5143015621
4_48   4_48   4   48 0.5230897401
4_49   4_49   4   49 0.5368007925
4_50   4_50   4   50 0.5661796148
4_51   4_51   4   51 0.5788620112
4_52   4_52   4   52 0.6230412417
4_53   4_53   4   53 0.6249633803
4_54   4_54   4   54 0.6311026101
4_55   4_55   4   55 0.6314469759
4_56   4_56   4   56 0.6610672071
4_57   4_57   4   57 0.6611086603
4_58   4_58   4   58 0.6618152977
4_59   4_59   4   59 0.6689019038
4_60   4_60   4   60 0.6713752543
4_61   4_61   4   61 0.6720634180
4_62   4_62   4   62 0.6768707512
4_63   4_63   4   63 0.6879729696
4_64   4_64   4   64 0.7027549669
4_65   4_65   4   65 0.7110308080
4_66   4_66   4   66 0.7190052902
4_67   4_67   4   67 0.7265234583
4_68   4_68   4   68 0.7393726817
4_69   4_69   4   69 0.7393829075
4_70   4_70   4   70 0.7428856302
4_71   4_71   4   71 0.7586378321
4_72   4_72   4   72 0.7636736985
4_73   4_73   4   73 0.7636847590
4_74   4_74   4   74 0.7659399379
4_75   4_75   4   75 0.7673790365
4_76   4_76   4   76 0.7711935610
4_77   4_77   4   77 0.7809714285
4_78   4_78   4   78 0.7865061170
4_79   4_79   4   79 0.7893353507
4_80   4_80   4   80 0.7918609544
4_81   4_81   4   81 0.8061675921
4_82   4_82   4   82 0.8104294839
4_83   4_83   4   83 0.8221287999
4_84   4_84   4   84 0.8244888510
4_85   4_85   4   85 0.8295543415
4_86   4_86   4   86 0.8498995148
4_87   4_87   4   87 0.8793448957
4_88   4_88   4   88 0.8839645072
4_89   4_89   4   89 0.9015264034
4_90   4_90   4   90 0.9049238717
4_91   4_91   4   91 0.9155315642
4_92   4_92   4   92 0.9162221332
4_93   4_93   4   93 0.9230823232
4_94   4_94   4   94 0.9266435664
4_95   4_95   4   95 0.9388967960
4_96   4_96   4   96 0.9410172980
4_97   4_97   4   97 0.9479380986
4_98   4_98   4   98 0.9497159280
4_99   4_99   4   99 0.9541394292
4_100 4_100   4  100 0.9779102918
5_1     5_1   5    1 0.0000000000
5_2     5_2   5    2 0.0060886337
5_3     5_3   5    3 0.0250593176
5_4     5_4   5    4 0.0271347532
5_5     5_5   5    5 0.0312867337
5_6     5_6   5    6 0.0539463117
5_7     5_7   5    7 0.0827549716
5_8     5_8   5    8 0.1002663257
5_9     5_9   5    9 0.1047184997
5_10   5_10   5   10 0.1105350739
5_11   5_11   5   11 0.1152704149
5_12   5_12   5   12 0.1379278101
5_13   5_13   5   13 0.1442130732
5_14   5_14   5   14 0.1713050204
5_15   5_15   5   15 0.2044335191
5_16   5_16   5   16 0.2135640830
5_17   5_17   5   17 0.2147022919
5_18   5_18   5   18 0.2155743155
5_19   5_19   5   19 0.2159894896
5_20   5_20   5   20 0.2276636883
5_21   5_21   5   21 0.2276823383
5_22   5_22   5   22 0.2412490482
5_23   5_23   5   23 0.2455072787
5_24   5_24   5   24 0.2568995213
5_25   5_25   5   25 0.2588479296
5_26   5_26   5   26 0.2615279635
5_27   5_27   5   27 0.2803512314
5_28   5_28   5   28 0.2828139633
5_29   5_29   5   29 0.2832239616
5_30   5_30   5   30 0.3006263503
5_31   5_31   5   31 0.3056797302
5_32   5_32   5   32 0.3104358661
5_33   5_33   5   33 0.3129433876
5_34   5_34   5   34 0.3216692723
5_35   5_35   5   35 0.3429150562
5_36   5_36   5   36 0.3780645881
5_37   5_37   5   37 0.3989878852
5_38   5_38   5   38 0.4021754865
5_39   5_39   5   39 0.4109402692
5_40   5_40   5   40 0.4154954046
5_41   5_41   5   41 0.4293330439
5_42   5_42   5   42 0.4448340287
5_43   5_43   5   43 0.4569708978
5_44   5_44   5   44 0.4615705587
5_45   5_45   5   45 0.4744353628
5_46   5_46   5   46 0.4811237506
5_47   5_47   5   47 0.4850809668
5_48   5_48   5   48 0.4986981213
5_49   5_49   5   49 0.5035714810
5_50   5_50   5   50 0.5169856444
5_51   5_51   5   51 0.5217335660
5_52   5_52   5   52 0.5255296888
5_53   5_53   5   53 0.5345028092
5_54   5_54   5   54 0.5364666912
5_55   5_55   5   55 0.5594113662
5_56   5_56   5   56 0.5751493556
5_57   5_57   5   57 0.5902678440
5_58   5_58   5   58 0.5991872913
5_59   5_59   5   59 0.6045948034
5_60   5_60   5   60 0.6157164701
5_61   5_61   5   61 0.6157552034
5_62   5_62   5   62 0.6188252125
5_63   5_63   5   63 0.6247577795
5_64   5_64   5   64 0.6280410151
5_65   5_65   5   65 0.6480095975
5_66   5_66   5   66 0.6563348959
5_67   5_67   5   67 0.6589672899
5_68   5_68   5   68 0.6674644980
5_69   5_69   5   69 0.6769275391
5_70   5_70   5   70 0.6788088523
5_71   5_71   5   71 0.6836269086
5_72   5_72   5   72 0.6887608158
5_73   5_73   5   73 0.7116735309
5_74   5_74   5   74 0.7173994846
5_75   5_75   5   75 0.7366725201
5_76   5_76   5   76 0.7379164125
5_77   5_77   5   77 0.7549740871
5_78   5_78   5   78 0.7662831048
5_79   5_79   5   79 0.7671912454
5_80   5_80   5   80 0.7806847844
5_81   5_81   5   81 0.7889908212
5_82   5_82   5   82 0.7986584099
5_83   5_83   5   83 0.8192949024
5_84   5_84   5   84 0.8226620115
5_85   5_85   5   85 0.8262408744
5_86   5_86   5   86 0.8262550312
5_87   5_87   5   87 0.8464299132
5_88   5_88   5   88 0.8499987317
5_89   5_89   5   89 0.8959359020
5_90   5_90   5   90 0.8959889912
5_91   5_91   5   91 0.9051739322
5_92   5_92   5   92 0.9081126947
5_93   5_93   5   93 0.9081754885
5_94   5_94   5   94 0.9144869797
5_95   5_95   5   95 0.9324247846
5_96   5_96   5   96 0.9550111631
5_97   5_97   5   97 0.9827952267
5_98   5_98   5   98 0.9853737082
5_99   5_99   5   99 0.9857209374
5_100 5_100   5  100 0.9907242892

######################################################################################################################################################
#Clean working environment
rmlist=ls()
#Set default pot layout
par(mfrow=c (1,1))
#load AlphasimR, simulate founder'S Genome
library(AlphaSimR)
founderGenome=quickHaplo(nInd = 100, nChr = 5, segSites = 100, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=100, mean=10, var=1)

#Base Population and their Phenotypes
basePop=newPop(founderGenome)
heritability=0.8
basePop=setPheno(basePop, H2=heritability)

#Let's summarize Phenotypes values in this base Population
hist(pheno(basePop))
meanP(basePop)
varP(basePop)

#Parent of Next Generation
basePopData=data.frame(id=basePop@id,Pheno=pheno(basePop),[, 1])
basePopData = data.frame(id = basePop@id, Pheno = pheno(basePop)[, 1])


#Save number of selected individuals for later use
nSelected=20

#Show top individuals by phenotype values
basePopData[order(basePopData$Pheno, decreasing = TRUE),][1: nSelected, ]

#Select Superior individuals
basePopSelected=selectInd(pop=basePop, nInd = nSelected, use="Pheno")

#Histogram of Phenotype values in Selected part of Base Population
hist(pheno(basePopSelected),xlim=range(pheno(basePop)), col="Purple")
#mean and variance of Selected Part of base Population
meanP(basePopSelected)
varP(basePopSelected)

#Cross Selected individuals
newPop=randCross(pop = basePopSelected, nCrosses = nInd(basePop))

#Phenotype the progeny
newPop=setPheno(newPop, h2=heritability)

#Histogram of Phenotype values
hist(pheno(newPop))

#mean and variance of phenotype values in new Population
meanP(newPop)
varP(newPop)

#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)
#Create population
pop = newPop(founderPop, simParam=SP)
> founderPop
An object of class "MapPop"
Ploidy: 2
Individuals: 10
Chromosomes: 1
Loci: 10
> class(founderPop)
[1] "MapPop"
attr(,"package")
[1] "AlphaSimR"
>
  > SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.8)
#Create population
pop = newPop(founderPop, simParam=SP)
> pop
An object of class "Pop"
Ploidy: 2
Individuals: 10
Chromosomes: 1
Loci: 10
Traits: 1
> class(pop)
[1] "Pop"
attr(,"package")
[1] "AlphaSimR"
> #Select top 5 (directional selection)
  pop2 = selectInd(pop, 5, simParam=SP)
> pop2
An object of class "Pop"
Ploidy: 2
Individuals: 5
Chromosomes: 1
Loci: 10
Traits: 1
> hist(pop@pheno); abline(v = pop@pheno, lwd = 2)
abline(v = pop2@pheno, col = "red", lwd = 2)
> hist(pop@pheno); abline(v = pop@pheno, col='red' lwd = 2)
abline(v = pop2@pheno, col = "blue", lwd = 2)
Error: unexpected symbol in " abline(v = pop@pheno, col='red' lwd"
> hist(pop@pheno); abline(v = pop@pheno, col='red' lwd = 2);abline(v = pop2@pheno, col = "blue", lwd = 2)
Error: unexpected symbol in " abline(v = pop@pheno, col='red' lwd"
> hist(pop@pheno); abline(v = pop@pheno, col='red', lwd = 2);abline(v = pop2@pheno, col = "blue", lwd = 2)
> pop@pheno
Trait1
[1,] -0.15961529
[2,]  1.30992835
[3,] -1.13572984
[4,]  1.27998485
[5,]  0.80345632
[6,] -1.21658966
[7,]  0.05757516
[8,] -1.56044263
[9,]  1.43295417
[10,]  0.81736719
> pop2@pheno
Trait1
[1,] 1.4329542
[2,] 1.3099284
[3,] 1.2799848
[4,] 0.8173672
[5,] 0.8034563
> plot(pop@pheno)
> boxplot(pop@pheno)
> mean(pop@pheno)
[1] 0.1628889
> sd(pop@pheno)
[1] 1.139672
> mean(pop2@pheno)
[1] 1.128738
> sor(pop@pheno)
Error in sor(pop@pheno) : could not find function "sor"
> sort(pop@pheno)
[1] -1.56044263 -1.21658966 -1.13572984 -0.15961529  0.05757516  0.80345632
[7]  0.81736719  1.27998485  1.30992835  1.43295417
> pop@
  pop@id       pop@mother   pop@sex      pop@gv       pop@ebv      pop@fixEff   pop@miscPop  pop@nChr     pop@nLoci
pop@iid      pop@father   pop@nTraits  pop@pheno    pop@gxe      pop@misc     pop@nInd     pop@ploidy   pop@geno
> pop@id
[1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10"
> pop@iid
[1]  1  2  3  4  5  6  7  8  9 10
> pop@mother
[1] "0" "0" "0" "0" "0" "0" "0" "0" "0" "0"
> pop@nTraits
[1] 1
> pop@gv
Trait1
[1,] -0.3281649
[2,]  0.8032636
[3,] -0.9902762
[4,]  0.7952615
[5,] -0.7244372
[6,] -1.3926573
[7,]  0.6814395
[8,] -1.2858176
[9,]  1.4101780
[10,]  1.0312106
> pop@ebv

[1,]
[2,]
[3,]
[4,]
[5,]
[6,]
[7,]
[8,]
[9,]
[10,]
> pop@geno
[[1]]
, , 1

[,1] [,2]
[1,]   7e   d8
[2,]   a5   d3

, , 2

[,1] [,2]
[1,]   3f   84
[2,]   0f   77

, , 3

[,1] [,2]
[1,]   e4   91
[2,]   05   87

, , 4

[,1] [,2]
[1,]   56   98
[2,]   21   3d

, , 5

[,1] [,2]
[1,]   d4   41
[2,]   a2   bc

, , 6

[,1] [,2]
[1,]   2d   e0
[2,]   94   6a

, , 7

[,1] [,2]
[1,]   81   8f
[2,]   bc   65

, , 8

[,1] [,2]
[1,]   f2   68
[2,]   d5   30

, , 9

[,1] [,2]
[1,]   ce   17
[2,]   03   1f

, , 10

[,1] [,2]
[1,]   3f   e7
[2,]   74   f2


> pop@gxe
[[1]]
NULL

> getPed(pop)
id mother father
1   1      0      0
2   2      0      0
3   3      0      0
4   4      0      0
5   5      0      0
6   6      0      0
7   7      0      0
8   8      0      0
9   9      0      0
10 10      0      0
> getQtlMap(
  trait=     sex=       simParam=
    > bv(pop)
  Trait1
  [1,] -0.3281649
  [2,]  0.8032636
  [3,] -0.9902762
  [4,]  0.7952615
  [5,] -0.7244372
  [6,] -1.3926573
  [7,]  0.6814395
  [8,] -1.2858176
  [9,]  1.4101780
  [10,]  1.0312106
  > bv(founderPop)
  Error in genParam(pop, simParam = simParam) :
    class(pop) == "Pop" is not TRUE
  > bv(pop)
  Trait1
  [1,] -0.3281649
  [2,]  0.8032636
  [3,] -0.9902762
  [4,]  0.7952615
  [5,] -0.7244372
  [6,] -1.3926573
  [7,]  0.6814395
  [8,] -1.2858176
  [9,]  1.4101780
  [10,]  1.0312106
  >   founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
  #Set simulation parameters
  SP = SimParam$new(founderPop)
  SP$addTraitA(10)
  > SP$setVarE(h2=1)
  > pop = newPop(founderPop, simParam=SP)
  > bv(pop)
  Trait1
  [1,] -1.56965534
  [2,] -0.45572314
  [3,]  0.74946995
  [4,] -1.12016274
  [5,] -0.42062259
  [6,]  1.92887549
  [7,] -0.48353159
  [8,]  0.05789664
  [9,]  1.16407768
  [10,]  0.14937563
  > pop@pheno
  Trait1
  [1,] -1.56965534
  [2,] -0.45572314
  [3,]  0.74946995
  [4,] -1.12016274
  [5,] -0.42062259
  [6,]  1.92887549
  [7,] -0.48353159
  [8,]  0.05789664
  [9,]  1.16407768
  [10,]  0.14937563
  
######################################################################################################################################################
  founderGenome=quickHaplo(nInd = 100, nChr = 5, segSites = 100, inbred=TRUE)
  SP=SimParam$new(founderGenome)
  SP$addTraitA(nQtlPerChr=100, mean=10, var=1)
  basePop=newPop(founderGenome)
  heritability=0.8
  basePop=setPheno(basePop, H2=heritability)
  hist(pheno(basePop))
  meanP(basePop)
  Trait1
  9.922858
  varP(basePop)
  Trait1
  Trait1 1.271647
  basePopData = data.frame(id = basePop@id, Pheno = pheno(basePop)[, 1])
  basePopData@
    +  basePopData@
  Error in basePopData@basePopData :
    no applicable method for `@` applied to an object of class "data.frame"
  basePop@id
  [1] "1"   "2"   "3"   "4"   "5"   "6"   "7"   "8"   "9"   "10"  "11"  "12"
  [13] "13"  "14"  "15"  "16"  "17"  "18"  "19"  "20"  "21"  "22"  "23"  "24"
  [25] "25"  "26"  "27"  "28"  "29"  "30"  "31"  "32"  "33"  "34"  "35"  "36"
  [37] "37"  "38"  "39"  "40"  "41"  "42"  "43"  "44"  "45"  "46"  "47"  "48"
  [49] "49"  "50"  "51"  "52"  "53"  "54"  "55"  "56"  "57"  "58"  "59"  "60"
  [61] "61"  "62"  "63"  "64"  "65"  "66"  "67"  "68"  "69"  "70"  "71"  "72"
  [73] "73"  "74"  "75"  "76"  "77"  "78"  "79"  "80"  "81"  "82"  "83"  "84"
  [85] "85"  "86"  "87"  "88"  "89"  "90"  "91"  "92"  "93"  "94"  "95"  "96"
  [97] "97"  "98"  "99"  "100"
  nSelected=20
  basePopSelected=selectInd(pop=basePop, nInd = nSelected, use="Pheno")
  basePopSelected@
    + basePop@pheno
  Error: no slot of name "basePop" for this object of class "Pop"
  basePopSelected@pheno
  Trait1
  [1,] 14.19049
  [2,] 12.16402
  [3,] 11.86842
  [4,] 11.83181
  [5,] 11.73691
  [6,] 11.72696
  [7,] 11.71569
  [8,] 11.70509
  [9,] 11.61273
  [10,] 11.36715
  [11,] 11.18546
  [12,] 11.17516
  [13,] 11.12368
  [14,] 11.12221
  [15,] 11.07114
  [16,] 11.06234
  [17,] 11.03134
  [18,] 11.02589
  [19,] 10.99818
  [20,] 10.93385
  > basePopData[order(basePopData$Pheno, decreasing = TRUE),][1: nSelected, ]
  id    Pheno
  78 78 14.19049
  33 33 12.16402
  22 22 11.86842
  45 45 11.83181
  86 86 11.73691
  19 19 11.72696
  46 46 11.71569
  82 82 11.70509
  37 37 11.61273
  13 13 11.36715
  17 17 11.18546
  2   2 11.17516
  47 47 11.12368
  56 56 11.12221
  58 58 11.07114
  69 69 11.06234
  65 65 11.03134
  36 36 11.02589
  28 28 10.99818
  75 75 10.93385
  > hist(pheno(basePopSelected),xlim=range(pheno(basePop)), col="Purple")
  > mean(basePopSelected)
  [1] NA
  Warning message:
    In mean.default(basePopSelected) :
    argument is not numeric or logical: returning NA
  > varP(basePopSelected
         + varP(basePopSelected
                Error: unexpected symbol in:
                  "varP(basePopSelected
varP"
                > varP(basePopSelected)
                Trait1
                Trait1 0.5006036
                > meanP(basePopSelected)
                Trait1
                11.53243
                > newPop=randCross(pop = basePopSelected, nCrosses = nInd(basePop))
                > newPop@
                  + newPop@id
                Error: no slot of name "newPop" for this object of class "Pop"
                > newPop=randCross(pop = basePopSelected, nCrosses = nInd(basePop))
                > newPop=setPheno(newPop, h2=heritability)
                > hist(pheno(newPop))
                > newPop@pheno
                Trait1
                [1,] 12.720072
                [2,] 12.195103
                [3,] 12.866510
                [4,] 12.818153
                [5,] 12.049293
                [6,] 13.332824
                [7,] 12.792024
                [8,] 12.048398
                [9,] 12.226003
                [10,] 11.889971
                [11,] 11.327076
                [12,] 11.464148
                [13,] 10.539877
                [14,] 11.335795
                [15,] 10.992399
                [16,] 10.854939
                [17,] 10.591051
                [18,] 10.915257
                [19,] 10.183128
                [20,] 10.421913
                [21,] 11.361788
                [22,] 10.855676
                [23,] 12.030421
                [24,] 10.756796
                [25,] 11.112474
                [26,] 10.544215
                [27,] 11.164580
                [28,] 12.156135
                [29,] 10.683520
                [30,] 11.035833
                [31,] 10.599834
                [32,] 10.186194
                [33,]  9.941344
                [34,]  9.839128
                [35,] 11.575421
                [36,] 11.639933
                [37,] 11.417394
                [38,] 11.284732
                [39,] 10.950791
                [40,] 10.366693
                [41,] 10.903755
                [42,]  9.843277
                [43,] 10.789706
                [44,] 10.475641
                [45,]  9.441800
                [46,] 10.958721
                [47,] 10.232643
                [48,] 11.144956
                [49,] 10.908165
                [50,]  9.612424
                [51,] 12.190430
                [52,] 11.587841
                [53,] 11.441056
                [54,] 11.781494
                [55,] 11.286567
                [56,] 12.129335
                [57,] 12.463551
                [58,] 11.871156
                [59,] 11.646805
                [60,] 12.398707
                [61,] 11.799549
                [62,] 10.768249
                [63,] 13.014060
                [64,] 10.361054
                [65,] 10.339241
                [66,] 11.584726
                [67,] 10.954228
                [68,] 10.305004
                [69,] 11.979703
                [70,] 11.390870
                [71,] 10.725582
                [72,] 10.833385
                [73,] 10.997601
                [74,] 10.368714
                [75,] 11.158332
                [76,] 11.091694
                [77,] 10.842622
                [78,] 10.586345
                [79,] 10.709314
                [80,] 10.464786
                [81,] 11.597700
                [82,] 11.007035
                [83,] 10.554276
                [84,] 11.767172
                [85,] 11.775283
                [86,] 10.692428
                [87,] 10.673805
                [88,] 11.151161
                [89,] 10.874428
                [90,] 10.880548
                [91,] 11.641017
                [92,] 10.961354
                [93,] 11.524702
                [94,] 10.844384
                [95,] 11.054144
                [96,] 10.462772
                [97,] 10.408731
                [98,] 11.727088
                [99,] 11.699972
                [100,]  9.855467
                > meanP(newPop)
                Trait1
                11.17603
                > varP(newPop)
                Trait1
                Trait1 0.6087644
plot(x=gv(basePop), y=pheno(basePop), xlab="Genetic value", ylab="Phenotype value")
cor(gv(basePop), pheno(basePop))
Trait1
Trait1 0.8838887

phenoRange=range(pheno(basePop, newPop))
                 
>
  
######################################################################################################################################################
founderGenome=quickHaplo(nInd = 10, nChr = 5, segSites = 100, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=10, mean=10, var=1)
basePop=newPop(founderGenome)

#Ext Haplotypes
basePopHaplo=pullSegSiteHaplo(basePop)
basePopHaplo [, 1:5]
basePopHaplo

#Ext Genotypes
basePopGeno=pullSegSiteGeno(basePop)
basePopGeno [, 1:5]

#To look variation and Resemblance between Genome of all individuals
-by simply counting no. of mutations each encode ancestral allele as 0 and mutation as 1.
rowSums(basePopHaplo)

#To count no. of mutations per Haplotype
rowSums(basePopGeno)

#Variation in progeny Genome
100 progeny for each crossprod(12, 13,32)

#count no. of mutations in parent genotypes
nPar=rowSums(basePopGeno)
nPar=nPar[1]
nPar=nPar[2]
nPar=nPar[3]

cross12=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))

#Ext Progeny genotypes

cross12Geno=pullSegSiteGeno(cross12)
#count no. of mutations in progeny genotypes
ncross12=rowSums(cross12Geno)
nPar=nPar[1]
nPar=nPar[2]
hist(c(cross12, cross13, cross23), xlim=rangeN, breaks=seq(from=rangeN[1], to=rangeN[2]),
     x(lab="Number of mutations"),
     abline(v=nPar1, col="blue", lwd=3),
     abline(v=nPar2, col="red", lwd=3),
     abline(v=nPar1+nPar2/2, col="black", lwd=3, lty=2)


cross13=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))

#Ext Progeny genotypes

cross13Geno=pullSegSiteGeno(cross13)
#count no. of mutations in progeny genotypes
ncross13=rowSums(cross13Geno)


cross23=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))


#Ext Progeny genotypes
#Ext Progeny genotypes

cross23Geno=pullSegSiteGeno(cross23)
pullSegSiteGeno(cross23)
#count no. of mutations in progeny genotypes
ncross23=rowSums(cross23Geno)
rowSums(cross23)


#All crosses together
hist(c(cross12, cross13, cross23), xlim=rangeN, breaks=seq(from=rangeN[1], to=rangeN[2]),
x(lab="Number of mutations"),
abline(v=nPar1, col="blue", lwd=3),
abline(v=nPar2, col="red", lwd=3),
abline(v=nPar3, col="green", lwd=3),
abline(v=nPar1+nPar2/2, col="black", lwd=3, lty=2)









heritability=0.2
basePop=setPheno(basePop, H2=heritability)
gv(basePop)
pheno(basePop)




###################################################################################################################################
founderGenome=quickHaplo(nInd = 10, nChr = 5, segSites = 100, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=10, mean=10, var=1)
basePop=newPop(founderGenome)
> basePopHaplo=pullSegSiteHaplo(basePop)
basePopHaplo [, 1:5]
basePopHaplo
1_1 1_2 1_3 1_4 1_5
1_1    0   0   1   0   1
1_2    0   0   1   0   1
2_1    1   0   1   1   1
2_2    1   0   1   1   1
3_1    0   1   0   1   0
3_2    0   1   0   1   0
4_1    1   1   0   1   1
4_2    1   1   0   1   1
5_1    0   1   1   1   0
5_2    0   1   1   1   0
6_1    1   1   1   0   1
6_2    1   1   1   0   1
7_1    1   0   1   0   0
7_2    1   0   1   0   0
8_1    0   0   0   1   0
8_2    0   0   0   1   0
9_1    0   1   1   0   1
9_2    0   1   1   0   1
10_1   0   1   0   0   0
10_2   0   1   0   0   0
1_1 1_2 1_3 1_4 1_5 1_6 1_7 1_8 1_9 1_10 1_11 1_12 1_13 1_14 1_15 1_16
1_1    0   0   1   0   1   0   0   1   0    0    1    0    1    1    0    0
1_2    0   0   1   0   1   0   0   1   0    0    1    0    1    1    0    0
2_1    1   0   1   1   1   1   0   0   0    1    1    1    1    0    0    1
2_2    1   0   1   1   1   1   0   0   0    1    1    1    1    0    0    1
3_1    0   1   0   1   0   1   1   0   0    0    1    1    0    0    1    0
3_2    0   1   0   1   0   1   1   0   0    0    1    1    0    0    1    0
4_1    1   1   0   1   1   0   1   1   0    1    0    1    1    1    1    0
4_2    1   1   0   1   1   0   1   1   0    1    0    1    1    1    1    0
5_1    0   1   1   1   0   1   1   1   0    0    0    1    0    1    0    1
5_2    0   1   1   1   0   1   1   1   0    0    0    1    0    1    0    1
6_1    1   1   1   0   1   0   1   1   0    1    0    0    0    1    0    1
6_2    1   1   1   0   1   0   1   1   0    1    0    0    0    1    0    1
7_1    1   0   1   0   0   0   0   1   0    1    0    1    0    0    0    1
7_2    1   0   1   0   0   0   0   1   0    1    0    1    0    0    0    1
8_1    0   0   0   1   0   1   1   0   0    0    0    1    1    0    1    0
8_2    0   0   0   1   0   1   1   0   0    0    0    1    1    0    1    0
9_1    0   1   1   0   1   1   0   1   0    1    1    1    1    1    1    0
9_2    0   1   1   0   1   1   0   1   0    1    1    1    1    1    1    0
10_1   0   1   0   0   0   1   1   1   1    1    0    1    1    0    0    0
10_2   0   1   0   0   0   1   1   1   1    1    0    1    1    0    0    0
1_17 1_18 1_19 1_20 1_21 1_22 1_23 1_24 1_25 1_26 1_27 1_28 1_29 1_30 1_31
1_1     1    0    1    1    1    0    0    1    0    1    0    0    0    1    0
1_2     1    0    1    1    1    0    0    1    0    1    0    0    0    1    0
2_1     1    1    0    1    1    0    1    0    1    0    1    1    1    1    0
2_2     1    1    0    1    1    0    1    0    1    0    1    1    1    1    0
3_1     1    0    1    0    0    1    0    1    0    1    0    0    1    0    1
3_2     1    0    1    0    0    1    0    1    0    1    0    0    1    0    1
4_1     1    0    1    0    0    0    1    0    1    0    0    0    1    1    0
4_2     1    0    1    0    0    0    1    0    1    0    0    0    1    1    0
5_1     0    1    0    1    1    1    0    0    0    0    1    0    0    0    1
5_2     0    1    0    1    1    1    0    0    0    0    1    0    0    0    1
6_1     0    1    1    1    0    0    1    0    1    0    1    0    1    1    1
6_2     0    1    1    1    0    0    1    0    1    0    1    0    1    1    1
7_1     0    0    1    0    1    1    0    1    1    0    1    1    0    1    1
7_2     0    0    1    0    1    1    0    1    1    0    1    1    0    1    1
8_1     1    0    1    0    1    1    0    1    1    1    0    1    1    1    0
8_2     1    0    1    0    1    1    0    1    1    1    0    1    1    1    0
9_1     1    0    0    1    1    0    0    1    1    0    0    0    1    1    1
9_2     1    0    0    1    1    0    0    1    1    0    0    0    1    1    1
10_1    0    1    0    1    1    1    1    1    1    1    1    1    1    1    1
10_2    0    1    0    1    1    1    1    1    1    1    1    1    1    1    1
1_32 1_33 1_34 1_35 1_36 1_37 1_38 1_39 1_40 1_41 1_42 1_43 1_44 1_45 1_46
1_1     0    1    1    1    0    1    0    1    1    1    1    0    0    1    0
1_2     0    1    1    1    0    1    0    1    1    1    1    0    0    1    0
2_1     1    0    1    1    1    0    0    0    0    0    0    1    1    1    1
2_2     1    0    1    1    1    0    0    0    0    0    0    1    1    1    1
3_1     1    0    1    1    1    0    0    1    1    0    1    0    1    1    0
3_2     1    0    1    1    1    0    0    1    1    0    1    0    1    1    0
4_1     1    0    1    0    0    1    1    0    1    1    1    0    1    0    1
4_2     1    0    1    0    0    1    1    0    1    1    1    0    1    0    1
5_1     0    1    0    1    0    1    0    1    1    0    0    0    1    0    1
5_2     0    1    0    1    0    1    0    1    1    0    0    0    1    0    1
6_1     1    1    1    1    1    0    0    1    1    1    0    1    0    1    0
6_2     1    1    1    1    1    0    0    1    1    1    0    1    0    1    0
7_1     0    0    0    0    1    0    0    0    0    0    1    0    1    0    1
7_2     0    0    0    0    1    0    0    0    0    0    1    0    1    0    1
8_1     0    0    0    1    0    1    0    0    1    1    0    0    0    1    1
8_2     0    0    0    1    0    1    0    0    1    1    0    0    0    1    1
9_1     1    1    1    1    1    0    1    0    0    1    1    1    0    0    1
9_2     1    1    1    1    1    0    1    0    0    1    1    1    0    0    1
10_1    1    0    0    1    1    1    0    0    1    1    1    1    1    1    1
10_2    1    0    0    1    1    1    0    0    1    1    1    1    1    1    1
1_47 1_48 1_49 1_50 1_51 1_52 1_53 1_54 1_55 1_56 1_57 1_58 1_59 1_60 1_61
1_1     1    0    0    0    0    0    1    1    1    0    1    1    1    0    0
1_2     1    0    0    0    0    0    1    1    1    0    1    1    1    0    0
2_1     0    1    0    0    0    0    0    0    0    1    0    1    0    0    1
2_2     0    1    0    0    0    0    0    0    0    1    0    1    0    0    1
3_1     0    1    1    1    0    1    1    0    1    0    0    1    0    0    1
3_2     0    1    1    1    0    1    1    0    1    0    0    1    0    0    1
4_1     1    1    1    0    0    0    1    1    1    0    1    0    0    0    0
4_2     1    1    1    0    0    0    1    1    1    0    1    0    0    0    0
5_1     1    0    0    1    1    0    1    1    1    1    1    0    0    1    0
5_2     1    0    0    1    1    0    1    1    1    1    1    0    0    1    0
6_1     1    0    1    0    0    0    1    1    0    1    0    1    0    1    0
6_2     1    0    1    0    0    0    1    1    0    1    0    1    0    1    0
7_1     0    0    1    0    1    0    0    0    1    0    1    0    1    0    0
7_2     0    0    1    0    1    0    0    0    1    0    1    0    1    0    0
8_1     0    1    1    1    0    1    1    0    0    0    1    0    0    1    0
8_2     0    1    1    1    0    1    1    0    0    0    1    0    0    1    0
9_1     1    0    0    1    1    0    0    0    1    1    0    1    1    0    0
9_2     1    0    0    1    1    0    0    0    1    1    0    1    1    0    0
10_1    0    0    1    1    1    0    1    1    1    0    1    0    1    1    0
10_2    0    0    1    1    1    0    1    1    1    0    1    0    1    1    0
1_62 1_63 1_64 1_65 1_66 1_67 1_68 1_69 1_70 1_71 1_72 1_73 1_74 1_75 1_76
1_1     1    0    1    1    1    0    0    0    0    0    0    1    0    0    1
1_2     1    0    1    1    1    0    0    0    0    0    0    1    0    0    1
2_1     0    0    0    1    1    1    1    0    1    1    1    0    1    1    0
2_2     0    0    0    1    1    1    1    0    1    1    1    0    1    1    0
3_1     1    0    0    1    0    0    1    0    1    0    0    0    0    0    0
3_2     1    0    0    1    0    0    1    0    1    0    0    0    0    0    0
4_1     0    1    1    1    1    1    0    1    1    0    0    1    1    0    0
4_2     0    1    1    1    1    1    0    1    1    0    0    1    1    0    0
5_1     1    0    1    0    1    1    0    1    0    1    0    0    0    0    1
5_2     1    0    1    0    1    1    0    1    0    1    0    0    0    0    1
6_1     0    1    0    1    1    0    1    1    1    1    0    0    0    1    1
6_2     0    1    0    1    1    0    1    1    1    1    0    0    0    1    1
7_1     0    1    0    1    1    1    0    1    0    1    0    1    0    1    1
7_2     0    1    0    1    1    1    0    1    0    1    0    1    0    1    1
8_1     1    0    1    1    0    1    1    0    1    0    0    0    1    0    1
8_2     1    0    1    1    0    1    1    0    1    0    0    0    1    0    1
9_1     1    1    0    0    0    0    0    1    0    0    0    1    0    1    0
9_2     1    1    0    0    0    0    0    1    0    0    0    1    0    1    0
10_1    1    0    1    0    1    1    0    0    1    0    0    1    0    1    1
10_2    1    0    1    0    1    1    0    0    1    0    0    1    0    1    1
1_77 1_78 1_79 1_80 1_81 1_82 1_83 1_84 1_85 1_86 1_87 1_88 1_89 1_90 1_91
1_1     1    1    1    1    1    1    0    0    1    0    0    0    1    1    1
1_2     1    1    1    1    1    1    0    0    1    0    0    0    1    1    1
2_1     1    0    0    1    0    1    0    1    0    0    0    0    1    1    0
2_2     1    0    0    1    0    1    0    1    0    0    0    0    1    1    0
3_1     0    1    0    0    0    0    0    1    0    1    1    1    1    0    0
3_2     0    1    0    0    0    0    0    1    0    1    1    1    1    0    0
4_1     0    1    1    0    0    1    0    1    1    0    1    1    1    0    1
4_2     0    1    1    0    0    1    0    1    1    0    1    1    1    0    1
5_1     1    0    1    1    0    1    1    1    1    1    0    1    0    0    1
5_2     1    0    1    1    0    1    1    1    1    1    0    1    0    0    1
6_1     0    1    1    1    1    1    1    0    0    0    1    0    0    0    0
6_2     0    1    1    1    1    1    1    0    0    0    1    0    0    0    0
7_1     0    0    1    0    0    1    0    1    0    1    1    1    0    0    1
7_2     0    0    1    0    0    1    0    1    0    1    1    1    0    0    1
8_1     0    1    0    0    0    1    1    1    0    0    0    0    1    0    0
8_2     0    1    0    0    0    1    1    1    0    0    0    0    1    0    0
9_1     0    1    0    1    1    0    1    1    1    0    0    0    1    1    1
9_2     0    1    0    1    1    0    1    1    1    0    0    0    1    1    1
10_1    0    1    1    0    0    1    1    1    1    1    1    0    0    1    1
10_2    0    1    1    0    0    1    1    1    1    1    1    0    0    1    1
1_92 1_93 1_94 1_95 1_96 1_97 1_98 1_99 1_100 2_1 2_2 2_3 2_4 2_5 2_6 2_7
1_1     1    1    0    0    0    0    0    1     1   1   1   0   1   1   0   0
1_2     1    1    0    0    0    0    0    1     1   1   1   0   1   1   0   0
2_1     0    1    0    1    1    0    0    0     1   0   1   0   1   1   0   0
2_2     0    1    0    1    1    0    0    0     1   0   1   0   1   1   0   0
3_1     1    1    1    0    0    1    1    0     0   0   0   1   1   1   1   1
3_2     1    1    1    0    0    1    1    0     0   0   0   1   1   1   1   1
4_1     0    1    1    1    1    0    1    0     0   0   1   1   1   0   0   1
4_2     0    1    1    1    1    0    1    0     0   0   1   1   1   0   0   1
5_1     1    0    0    1    0    0    1    1     0   1   0   0   0   1   0   0
5_2     1    0    0    1    0    0    1    1     0   1   0   0   0   1   0   0
6_1     0    1    0    0    0    0    0    1     1   1   0   0   0   1   1   0
6_2     0    1    0    0    0    0    0    1     1   1   0   0   0   1   1   0
7_1     1    1    1    0    1    1    0    1     1   0   0   1   1   1   1   0
7_2     1    1    1    0    1    1    0    1     1   0   0   1   1   1   1   0
8_1     1    1    1    1    1    0    1    1     0   1   1   0   0   1   0   0
8_2     1    1    1    1    1    0    1    1     0   1   1   0   0   1   0   0
9_1     1    1    0    0    1    0    0    0     1   0   0   0   0   1   1   0
9_2     1    1    0    0    1    0    0    0     1   0   0   0   0   1   1   0
10_1    1    1    1    0    1    0    0    1     0   1   0   1   1   0   1   1
10_2    1    1    1    0    1    0    0    1     0   1   0   1   1   0   1   1
2_8 2_9 2_10 2_11 2_12 2_13 2_14 2_15 2_16 2_17 2_18 2_19 2_20 2_21 2_22
1_1    1   0    0    0    0    0    0    0    1    1    1    1    0    1    0
1_2    1   0    0    0    0    0    0    0    1    1    1    1    0    1    0
2_1    1   0    1    0    1    0    1    0    0    0    1    0    1    0    1
2_2    1   0    1    0    1    0    1    0    0    0    1    0    1    0    1
3_1    0   1    1    1    0    1    0    0    0    0    1    0    0    0    0
3_2    0   1    1    1    0    1    0    0    0    0    1    0    0    0    0
4_1    0   0    1    0    0    0    1    0    0    0    0    0    0    1    0
4_2    0   0    1    0    0    0    1    0    0    0    0    0    0    1    0
5_1    0   0    0    1    0    0    0    0    0    1    1    0    1    0    0
5_2    0   0    0    1    0    0    0    0    0    1    1    0    1    0    0
6_1    1   0    1    0    1    0    1    0    0    0    1    1    1    1    0
6_2    1   0    1    0    1    0    1    0    0    0    1    1    1    1    0
7_1    1   0    1    1    0    1    1    0    1    1    0    1    0    1    0
7_2    1   0    1    1    0    1    1    0    1    1    0    1    0    1    0
8_1    1   0    1    1    0    1    0    0    1    1    0    0    1    1    1
8_2    1   0    1    1    0    1    0    0    1    1    0    0    1    1    1
9_1    1   1    1    0    1    0    0    1    1    1    0    1    0    1    1
9_2    1   1    1    0    1    0    0    1    1    1    0    1    0    1    1
10_1   0   0    0    1    0    0    1    1    0    1    0    0    1    1    1
10_2   0   0    0    1    0    0    1    1    0    1    0    0    1    1    1
2_23 2_24 2_25 2_26 2_27 2_28 2_29 2_30 2_31 2_32 2_33 2_34 2_35 2_36 2_37
1_1     0    0    1    1    1    0    1    1    1    0    0    1    0    1    0
1_2     0    0    1    1    1    0    1    1    1    0    0    1    0    1    0
2_1     0    0    0    1    1    0    1    1    0    1    1    0    1    0    0
2_2     0    0    0    1    1    0    1    1    0    1    1    0    1    0    0
3_1     1    1    0    1    0    0    1    0    1    0    1    0    1    0    1
3_2     1    1    0    1    0    0    1    0    1    0    1    0    1    0    1
4_1     0    0    1    1    0    1    0    0    1    0    1    1    0    0    0
4_2     0    0    1    1    0    1    0    0    1    0    1    1    0    0    0
5_1     0    1    1    1    1    1    1    0    1    1    0    1    1    0    1
5_2     0    1    1    1    1    1    1    0    1    1    0    1    1    0    1
6_1     1    1    0    1    0    1    1    0    1    0    1    0    1    0    0
6_2     1    1    0    1    0    1    1    0    1    0    1    0    1    0    0
7_1     0    0    0    0    1    0    0    1    1    0    0    0    1    0    1
7_2     0    0    0    0    1    0    0    1    1    0    0    0    1    0    1
8_1     1    0    0    1    1    0    1    1    1    0    0    0    0    1    1
8_2     1    0    0    1    1    0    1    1    1    0    0    0    0    1    1
9_1     1    1    0    1    0    0    1    0    1    0    1    1    0    1    0
9_2     1    1    0    1    0    0    1    0    1    0    1    1    0    1    0
10_1    1    1    0    0    0    1    0    1    1    1    1    0    0    0    0
10_2    1    1    0    0    0    1    0    1    1    1    1    0    0    0    0
2_38 2_39 2_40 2_41 2_42 2_43 2_44 2_45 2_46 2_47 2_48 2_49 2_50 2_51 2_52
1_1     0    0    0    0    1    0    0    1    1    1    0    0    0    1    1
1_2     0    0    0    0    1    0    0    1    1    1    0    0    0    1    1
2_1     1    1    1    1    1    1    0    1    0    0    1    1    0    0    1
2_2     1    1    1    1    1    1    0    1    0    0    1    1    0    0    1
3_1     0    0    1    0    1    0    1    0    0    0    0    1    0    1    1
3_2     0    0    1    0    1    0    1    0    0    0    0    1    0    1    1
4_1     0    1    1    0    0    1    0    1    1    0    0    1    0    1    0
4_2     0    1    1    0    0    1    0    1    1    0    0    1    0    1    0
5_1     1    0    0    1    0    1    0    1    1    1    0    1    1    0    0
5_2     1    0    0    1    0    1    0    1    1    1    0    1    1    0    0
6_1     0    0    1    0    0    0    0    1    1    1    1    0    0    1    1
6_2     0    0    1    0    0    0    0    1    1    1    1    0    0    1    1
7_1     1    1    0    1    1    0    1    0    0    0    1    0    1    1    0
7_2     1    1    0    1    1    0    1    0    0    0    1    0    1    1    0
8_1     0    1    0    0    1    0    1    1    0    1    0    0    1    0    1
8_2     0    1    0    0    1    0    1    1    0    1    0    0    1    0    1
9_1     1    0    1    0    0    1    0    0    1    1    1    1    1    1    0
9_2     1    0    1    0    0    1    0    0    1    1    1    1    1    1    0
10_1    1    1    0    1    1    0    0    1    0    0    0    1    0    1    0
10_2    1    1    0    1    1    0    0    1    0    0    0    1    0    1    0
2_53 2_54 2_55 2_56 2_57 2_58 2_59 2_60 2_61 2_62 2_63 2_64 2_65 2_66 2_67
1_1     1    1    0    1    0    1    1    0    0    1    1    0    0    0    1
1_2     1    1    0    1    0    1    1    0    0    1    1    0    0    0    1
2_1     0    1    1    1    1    0    0    0    0    1    0    0    1    0    1
2_2     0    1    1    1    1    0    0    0    0    1    0    0    1    0    1
3_1     0    0    0    1    0    0    1    0    0    0    1    1    1    1    1
3_2     0    0    0    1    0    0    1    0    0    0    1    1    1    1    1
4_1     0    0    1    1    1    0    1    0    0    1    0    1    0    1    1
4_2     0    0    1    1    1    0    1    0    0    1    0    1    0    1    1
5_1     1    1    1    1    1    1    0    1    1    0    0    1    1    0    1
5_2     1    1    1    1    1    1    0    1    1    0    0    1    1    0    1
6_1     0    0    0    1    1    0    0    1    0    1    0    1    0    0    0
6_2     0    0    0    1    1    0    0    1    0    1    0    1    0    0    0
7_1     0    0    1    0    0    0    1    1    1    1    1    1    0    1    0
7_2     0    0    1    0    0    0    1    1    1    1    1    1    0    1    0
8_1     1    1    1    1    0    0    1    1    1    0    0    0    1    1    0
8_2     1    1    1    1    0    0    1    1    1    0    0    0    1    1    0
9_1     1    0    0    0    1    1    0    1    0    1    0    0    0    1    0
9_2     1    0    0    0    1    1    0    1    0    1    0    0    0    1    0
10_1    1    0    0    1    0    0    1    1    1    1    1    1    0    1    0
10_2    1    0    0    1    0    0    1    1    1    1    1    1    0    1    0
2_68 2_69 2_70 2_71 2_72 2_73 2_74 2_75 2_76 2_77 2_78 2_79 2_80 2_81 2_82
1_1     1    1    0    0    0    0    0    1    0    0    0    1    0    1    1
1_2     1    1    0    0    0    0    0    1    0    0    0    1    0    1    1
2_1     0    0    0    0    1    0    1    0    0    1    0    1    0    1    0
2_2     0    0    0    0    1    0    1    0    0    1    0    1    0    1    0
3_1     1    0    0    1    1    0    1    0    0    1    1    0    1    0    1
3_2     1    0    0    1    1    0    1    0    0    1    1    0    1    0    1
4_1     1    1    1    1    0    1    1    0    0    1    0    0    0    0    0
4_2     1    1    1    1    0    1    1    0    0    1    0    0    0    0    0
5_1     1    0    0    0    1    0    1    1    0    1    1    0    0    0    1
5_2     1    0    0    0    1    0    1    1    0    1    1    0    0    0    1
6_1     0    0    0    1    1    0    1    0    1    1    0    1    0    1    0
6_2     0    0    0    1    1    0    1    0    1    1    0    1    0    1    0
7_1     1    1    0    0    1    1    0    1    1    1    0    1    1    0    0
7_2     1    1    0    0    1    1    0    1    1    1    0    1    1    0    0
8_1     1    1    1    0    0    0    0    1    0    1    1    1    1    0    0
8_2     1    1    1    0    0    0    0    1    0    1    1    1    1    0    0
9_1     1    0    1    0    1    1    1    0    0    0    0    0    0    0    1
9_2     1    0    1    0    1    1    1    0    0    0    0    0    0    0    1
10_1    0    1    0    1    0    0    0    1    1    0    0    1    0    0    0
10_2    0    1    0    1    0    0    0    1    1    0    0    1    0    0    0
2_83 2_84 2_85 2_86 2_87 2_88 2_89 2_90 2_91 2_92 2_93 2_94 2_95 2_96 2_97
1_1     1    1    1    1    1    0    1    0    1    1    1    0    0    1    0
1_2     1    1    1    1    1    0    1    0    1    1    1    0    0    1    0
2_1     0    0    1    1    0    0    0    1    0    1    1    1    1    1    0
2_2     0    0    1    1    0    0    0    1    0    1    1    1    1    1    0
3_1     0    1    0    0    0    0    0    1    0    0    0    0    1    1    0
3_2     0    1    0    0    0    0    0    1    0    0    0    0    1    1    0
4_1     0    0    0    1    1    0    0    0    0    1    0    1    1    0    0
4_2     0    0    0    1    1    0    0    0    0    1    0    1    1    0    0
5_1     0    1    0    1    0    1    1    0    0    1    0    1    1    1    0
5_2     0    1    0    1    0    1    1    0    0    1    0    1    1    1    0
6_1     0    1    0    1    0    1    0    1    1    1    0    1    0    1    1
6_2     0    1    0    1    0    1    0    1    1    1    0    1    0    1    1
7_1     1    0    0    0    0    1    1    1    1    1    0    1    0    1    0
7_2     1    0    0    0    0    1    1    1    1    1    0    1    0    1    0
8_1     0    0    1    1    1    0    1    0    0    0    1    1    0    1    0
8_2     0    0    1    1    1    0    1    0    0    0    1    1    0    1    0
9_1     0    0    0    0    1    1    1    1    0    0    0    0    0    0    0
9_2     0    0    0    0    1    1    1    1    0    0    0    0    0    0    0
10_1    1    1    1    1    1    1    0    1    0    0    0    0    0    0    0
10_2    1    1    1    1    1    1    0    1    0    0    0    0    0    0    0
2_98 2_99 2_100 3_1 3_2 3_3 3_4 3_5 3_6 3_7 3_8 3_9 3_10 3_11 3_12 3_13
1_1     0    1     1   0   0   1   0   1   1   0   1   1    1    1    1    0
1_2     0    1     1   0   0   1   0   1   1   0   1   1    1    1    1    0
2_1     1    0     0   0   0   0   0   1   0   1   1   0    0    0    1    0
2_2     1    0     0   0   0   0   0   1   0   1   1   0    0    0    1    0
3_1     1    1     0   1   0   1   0   0   0   1   1   1    1    1    0    0
3_2     1    1     0   1   0   1   0   0   0   1   1   1    1    1    0    0
4_1     1    0     1   0   0   0   1   1   1   0   1   1    0    1    1    0
4_2     1    0     1   0   0   0   1   1   1   0   1   1    0    1    1    0
5_1     1    0     1   0   0   0   0   0   0   1   1   1    1    1    1    0
5_2     1    0     1   0   0   0   0   0   0   1   1   1    1    1    1    0
6_1     1    0     0   1   1   0   0   1   1   1   0   0    0    1    1    0
6_2     1    0     0   1   1   0   0   1   1   1   0   0    0    1    1    0
7_1     1    1     0   0   1   0   1   0   1   1   1   0    0    1    0    1
7_2     1    1     0   0   1   0   1   0   1   1   1   0    0    1    0    1
8_1     1    1     0   0   0   0   1   1   1   1   0   0    0    0    0    1
8_2     1    1     0   0   0   0   1   1   1   1   0   0    0    0    0    1
9_1     1    1     1   1   0   0   0   1   1   1   0   0    1    0    0    0
9_2     1    1     1   1   0   0   0   1   1   1   0   0    1    0    0    0
10_1    1    0     1   1   1   1   0   0   1   0   1   1    0    0    0    1
10_2    1    0     1   1   1   1   0   0   1   0   1   1    0    0    0    1
3_14 3_15 3_16 3_17 3_18 3_19 3_20 3_21 3_22 3_23 3_24 3_25 3_26 3_27 3_28
1_1     1    1    1    0    1    0    1    0    0    0    1    1    0    1    1
1_2     1    1    1    0    1    0    1    0    0    0    1    1    0    1    1
2_1     0    0    0    1    0    1    1    1    1    1    0    0    0    0    0
2_2     0    0    0    1    0    1    1    1    1    1    0    0    0    0    0
3_1     1    1    0    1    1    1    1    0    0    0    1    1    1    1    0
3_2     1    1    0    1    1    1    1    0    0    0    1    1    1    1    0
4_1     1    1    0    1    0    1    1    1    0    0    0    1    1    0    0
4_2     1    1    0    1    0    1    1    1    0    0    0    1    1    0    0
5_1     0    1    1    0    0    1    0    1    0    0    0    0    0    1    1
5_2     0    1    1    0    0    1    0    1    0    0    0    0    0    1    1
6_1     1    0    1    0    0    0    0    1    1    1    0    0    1    0    1
6_2     1    0    1    0    0    0    0    1    1    1    0    0    1    0    1
7_1     1    0    0    1    0    0    1    0    0    0    1    0    0    1    0
7_2     1    0    0    1    0    0    1    0    0    0    1    0    0    1    0
8_1     1    1    0    1    1    0    1    1    1    1    1    1    1    0    1
8_2     1    1    0    1    1    0    1    1    1    1    1    1    1    0    1
9_1     0    1    0    0    1    0    0    0    0    1    0    0    1    1    0
9_2     0    1    0    0    1    0    0    0    0    1    0    0    1    1    0
10_1    1    1    1    1    1    0    0    0    0    1    0    0    1    0    1
10_2    1    1    1    1    1    0    0    0    0    1    0    0    1    0    1
3_29 3_30 3_31 3_32 3_33 3_34 3_35 3_36 3_37 3_38 3_39 3_40 3_41 3_42 3_43
1_1     1    0    1    1    1    0    1    1    0    1    0    0    0    1    0
1_2     1    0    1    1    1    0    1    1    0    1    0    0    0    1    0
2_1     0    1    0    1    1    1    1    1    1    1    1    0    1    0    1
2_2     0    1    0    1    1    1    1    1    1    1    1    0    1    0    1
3_1     1    1    0    1    1    1    0    1    0    0    0    0    0    1    1
3_2     1    1    0    1    1    1    0    1    0    0    0    0    0    1    1
4_1     1    0    0    1    1    1    0    1    0    1    0    1    0    0    1
4_2     1    0    0    1    1    1    0    1    0    1    0    1    0    0    1
5_1     0    0    0    0    0    1    1    0    1    1    1    1    0    0    1
5_2     0    0    0    0    0    1    1    0    1    1    1    1    0    0    1
6_1     1    0    1    1    1    0    1    0    0    1    1    0    0    0    1
6_2     1    0    1    1    1    0    1    0    0    1    1    0    0    0    1
7_1     0    1    0    1    0    0    1    0    1    1    0    1    0    1    1
7_2     0    1    0    1    0    0    1    0    1    1    0    1    0    1    1
8_1     0    1    1    0    0    1    1    0    0    1    0    1    0    1    0
8_2     0    1    1    0    0    1    1    0    0    1    0    1    0    1    0
9_1     0    1    1    0    1    0    0    0    0    1    0    1    0    1    0
9_2     0    1    1    0    1    0    0    0    0    1    0    1    0    1    0
10_1    0    1    0    1    0    1    0    0    0    0    1    0    1    1    1
10_2    0    1    0    1    0    1    0    0    0    0    1    0    1    1    1
3_44 3_45 3_46 3_47 3_48 3_49 3_50 3_51 3_52 3_53 3_54 3_55 3_56 3_57 3_58
1_1     0    1    0    0    0    0    0    0    1    1    0    1    0    1    1
1_2     0    1    0    0    0    0    0    0    1    1    0    1    0    1    1
2_1     1    1    0    1    0    1    0    0    0    1    1    0    0    0    1
2_2     1    1    0    1    0    1    0    0    0    1    1    0    0    0    1
3_1     0    1    0    1    0    1    0    0    0    1    0    1    1    0    1
3_2     0    1    0    1    0    1    0    0    0    1    0    1    1    0    1
4_1     1    0    1    1    0    0    1    1    0    0    0    1    1    0    1
4_2     1    0    1    1    0    0    1    1    0    0    0    1    1    0    1
5_1     1    1    0    1    1    1    0    0    0    1    1    0    1    0    1
5_2     1    1    0    1    1    1    0    0    0    1    1    0    1    0    1
6_1     1    0    1    1    0    0    1    1    0    1    0    0    1    0    0
6_2     1    0    1    1    0    0    1    1    0    1    0    0    1    0    0
7_1     0    1    1    0    0    0    0    0    1    0    0    0    0    1    0
7_2     0    1    1    0    0    0    0    0    1    0    0    0    0    1    0
8_1     1    1    1    0    1    1    1    0    0    1    1    0    0    0    1
8_2     1    1    1    0    1    1    1    0    0    1    1    0    0    0    1
9_1     0    1    0    0    1    0    0    1    1    0    1    0    1    0    1
9_2     0    1    0    0    1    0    0    1    1    0    1    0    1    0    1
10_1    1    0    1    1    0    0    0    0    1    1    1    1    1    1    1
10_2    1    0    1    1    0    0    0    0    1    1    1    1    1    1    1
3_59 3_60 3_61 3_62 3_63 3_64 3_65 3_66 3_67 3_68 3_69 3_70 3_71 3_72 3_73
1_1     1    1    1    1    0    0    0    1    1    0    1    1    0    1    0
1_2     1    1    1    1    0    0    0    1    1    0    1    1    0    1    0
2_1     1    1    0    0    1    0    0    0    0    0    1    1    1    0    0
2_2     1    1    0    0    1    0    0    0    0    0    1    1    1    0    0
3_1     1    1    1    0    1    0    0    1    1    1    1    1    1    0    1
3_2     1    1    1    0    1    0    0    1    1    1    1    1    1    0    1
4_1     1    1    0    1    1    1    1    1    1    1    0    0    1    1    1
4_2     1    1    0    1    1    1    1    1    1    1    0    0    1    1    1
5_1     0    1    1    0    0    0    0    1    1    1    0    1    1    1    0
5_2     0    1    1    0    0    0    0    1    1    1    0    1    1    1    0
6_1     1    0    1    1    0    0    0    1    0    1    1    1    1    0    1
6_2     1    0    1    1    0    0    0    1    0    1    1    1    1    0    1
7_1     1    0    0    1    1    1    1    0    1    1    0    1    1    1    0
7_2     1    0    0    1    1    1    1    0    1    1    0    1    1    1    0
8_1     1    1    0    1    0    0    1    0    1    0    0    1    0    1    1
8_2     1    1    0    1    0    0    1    0    1    0    0    1    0    1    1
9_1     0    1    0    1    0    1    1    0    0    1    1    0    0    0    0
9_2     0    1    0    1    0    1    1    0    0    1    1    0    0    0    0
10_1    0    0    0    1    1    0    1    1    1    1    0    0    0    0    1
10_2    0    0    0    1    1    0    1    1    1    1    0    0    0    0    1
3_74 3_75 3_76 3_77 3_78 3_79 3_80 3_81 3_82 3_83 3_84 3_85 3_86 3_87 3_88
1_1     1    0    1    1    1    1    0    0    1    1    0    0    0    0    1
1_2     1    0    1    1    1    1    0    0    1    1    0    0    0    0    1
2_1     0    1    0    1    1    1    0    0    1    1    1    1    1    1    1
2_2     0    1    0    1    1    1    0    0    1    1    1    1    1    1    1
3_1     0    0    0    0    0    0    0    0    1    1    0    0    0    0    1
3_2     0    0    0    0    0    0    0    0    1    1    0    0    0    0    1
4_1     1    0    0    1    1    0    1    0    0    1    0    1    0    0    1
4_2     1    0    0    1    1    0    1    0    0    1    0    1    0    0    1
5_1     1    0    1    0    1    0    1    1    1    1    0    1    1    0    1
5_2     1    0    1    0    1    0    1    1    1    1    0    1    1    0    1
6_1     1    1    0    0    1    0    1    0    0    0    1    1    0    0    0
6_2     1    1    0    0    1    0    1    0    0    0    1    1    0    0    0
7_1     1    1    0    1    0    1    1    1    0    0    0    0    1    1    1
7_2     1    1    0    1    0    1    1    1    0    0    0    0    1    1    1
8_1     0    1    0    0    1    1    0    1    0    0    1    1    0    0    0
8_2     0    1    0    0    1    1    0    1    0    0    1    1    0    0    0
9_1     1    0    1    0    1    0    0    1    1    1    1    1    1    0    1
9_2     1    0    1    0    1    0    0    1    1    1    1    1    1    0    1
10_1    0    0    1    1    1    0    0    1    0    1    1    1    1    0    1
10_2    0    0    1    1    1    0    0    1    0    1    1    1    1    0    1
3_89 3_90 3_91 3_92 3_93 3_94 3_95 3_96 3_97 3_98 3_99 3_100 4_1 4_2 4_3
1_1     1    0    1    1    0    0    1    0    1    1    0     1   0   1   0
1_2     1    0    1    1    0    0    1    0    1    1    0     1   0   1   0
2_1     0    1    0    0    1    1    1    1    0    1    1     1   1   1   1
2_2     0    1    0    0    1    1    1    1    0    1    1     1   1   1   1
3_1     0    1    1    1    1    1    1    0    1    1    1     1   1   0   0
3_2     0    1    1    1    1    1    1    0    1    1    1     1   1   0   0
4_1     0    1    0    1    0    1    0    1    1    0    0     0   0   0   1
4_2     0    1    0    1    0    1    0    1    1    0    0     0   0   0   1
5_1     0    0    1    1    1    0    1    0    0    1    1     1   0   0   0
5_2     0    0    1    1    1    0    1    0    0    1    1     1   0   0   0
6_1     0    1    1    1    1    1    0    0    1    1    1     1   1   0   1
6_2     0    1    1    1    1    1    0    0    1    1    1     1   1   0   1
7_1     1    0    1    1    1    1    0    1    0    1    0     1   0   1   0
7_2     1    0    1    1    1    1    0    1    0    1    0     1   0   1   0
8_1     0    0    1    1    0    1    0    0    0    0    0     1   0   1   1
8_2     0    0    1    1    0    1    0    0    0    0    0     1   0   1   1
9_1     1    0    0    0    0    0    0    1    1    1    1     1   0   1   0
9_2     1    0    0    0    0    0    0    1    1    1    1     1   0   1   0
10_1    1    0    0    1    0    0    1    1    0    0    0     0   1   0   1
10_2    1    0    0    1    0    0    1    1    0    0    0     0   1   0   1
4_4 4_5 4_6 4_7 4_8 4_9 4_10 4_11 4_12 4_13 4_14 4_15 4_16 4_17 4_18 4_19
1_1    0   0   1   0   1   1    0    0    1    0    0    0    0    0    1    0
1_2    0   0   1   0   1   1    0    0    1    0    0    0    0    0    1    0
2_1    0   1   0   0   1   0    0    0    0    1    0    0    1    0    0    1
2_2    0   1   0   0   1   0    0    0    0    1    0    0    1    0    0    1
3_1    0   0   1   0   1   0    1    1    0    1    0    0    0    0    1    0
3_2    0   0   1   0   1   0    1    1    0    1    0    0    0    0    1    0
4_1    1   0   0   1   1   1    1    1    0    1    0    0    0    0    1    0
4_2    1   0   0   1   1   1    1    1    0    1    0    0    0    0    1    0
5_1    0   0   0   0   0   0    0    1    1    0    1    1    0    0    0    1
5_2    0   0   0   0   0   0    0    1    1    0    1    1    0    0    0    1
6_1    0   0   0   1   1   1    1    0    0    0    0    1    1    1    0    1
6_2    0   0   0   1   1   1    1    0    0    0    0    1    1    1    0    1
7_1    0   0   1   1   0   0    1    0    1    0    0    1    0    1    1    1
7_2    0   0   1   1   0   0    1    0    1    0    0    1    0    1    1    1
8_1    1   0   0   0   1   1    0    0    1    1    0    1    1    0    0    0
8_2    1   0   0   0   1   1    0    0    1    1    0    1    1    0    0    0
9_1    0   0   1   0   1   0    1    0    0    1    1    1    1    1    0    0
9_2    0   0   1   0   1   0    1    0    0    1    1    1    1    1    0    0
10_1   1   0   0   0   1   1    1    1    0    1    1    0    1    1    1    1
10_2   1   0   0   0   1   1    1    1    0    1    1    0    1    1    1    1
4_20 4_21 4_22 4_23 4_24 4_25 4_26 4_27 4_28 4_29 4_30 4_31 4_32 4_33 4_34
1_1     1    0    1    0    0    1    0    1    0    0    1    0    1    0    1
1_2     1    0    1    0    0    1    0    1    0    0    1    0    1    0    1
2_1     1    1    1    1    0    0    0    0    0    0    0    0    1    1    1
2_2     1    1    1    1    0    0    0    0    0    0    0    0    1    1    1
3_1     0    0    1    1    0    0    0    0    0    0    0    1    0    1    0
3_2     0    0    1    1    0    0    0    0    0    0    0    1    0    1    0
4_1     1    0    0    0    0    1    1    0    0    0    1    0    0    1    0
4_2     1    0    0    0    0    1    1    0    0    0    1    0    0    1    0
5_1     1    1    0    0    0    0    0    0    0    0    1    0    0    0    1
5_2     1    1    0    0    0    0    0    0    0    0    1    0    0    0    1
6_1     0    0    0    1    0    1    1    0    0    1    1    1    1    0    1
6_2     0    0    0    1    0    1    1    0    0    1    1    1    1    0    1
7_1     0    0    0    0    0    1    1    1    1    1    0    1    0    1    0
7_2     0    0    0    0    0    1    1    1    1    1    0    1    0    1    0
8_1     1    1    1    0    1    0    1    1    1    1    1    0    1    0    0
8_2     1    1    1    0    1    0    1    1    1    1    1    0    1    0    0
9_1     0    0    0    0    0    1    1    0    0    1    1    1    0    1    0
9_2     0    0    0    0    0    1    1    0    0    1    1    1    0    1    0
10_1    0    1    1    1    0    1    1    0    1    0    1    0    1    1    1
10_2    0    1    1    1    0    1    1    0    1    0    1    0    1    1    1
4_35 4_36 4_37 4_38 4_39 4_40 4_41 4_42 4_43 4_44 4_45 4_46 4_47 4_48 4_49
1_1     0    1    0    1    0    0    1    1    1    0    1    1    1    1    0
1_2     0    1    0    1    0    0    1    1    1    0    1    1    1    1    0
2_1     1    1    0    0    1    0    0    0    1    1    0    1    1    0    1
2_2     1    1    0    0    1    0    0    0    1    1    0    1    1    0    1
3_1     1    0    0    0    0    1    0    1    1    1    1    1    1    0    1
3_2     1    0    0    0    0    1    0    1    1    1    1    1    1    0    1
4_1     1    0    0    0    1    0    1    0    1    1    0    1    0    1    1
4_2     1    0    0    0    1    0    1    0    1    1    0    1    0    1    1
5_1     0    1    1    0    1    0    1    1    1    1    1    0    1    0    1
5_2     0    1    1    0    1    0    1    1    1    1    1    0    1    0    1
6_1     1    0    0    1    1    0    0    0    0    0    1    0    1    1    1
6_2     1    0    0    1    1    0    0    0    0    0    1    0    1    1    1
7_1     0    1    1    1    1    0    0    1    0    1    1    0    1    0    1
7_2     0    1    1    1    1    0    0    1    0    1    1    0    1    0    1
8_1     1    1    1    0    1    1    1    1    1    1    1    0    1    0    1
8_2     1    1    1    0    1    1    1    1    1    1    1    0    1    0    1
9_1     0    1    0    0    0    0    1    0    1    1    1    1    1    1    1
9_2     0    1    0    0    0    0    1    0    1    1    1    1    1    1    1
10_1    1    1    0    1    0    0    0    0    0    1    1    0    0    0    1
10_2    1    1    0    1    0    0    0    0    0    1    1    0    0    0    1
4_50 4_51 4_52 4_53 4_54 4_55 4_56 4_57 4_58 4_59 4_60 4_61 4_62 4_63 4_64
1_1     1    0    1    0    0    1    1    0    1    1    1    1    0    1    1
1_2     1    0    1    0    0    1    1    0    1    1    1    1    0    1    1
2_1     1    1    0    1    1    0    0    0    1    0    1    0    1    0    0
2_2     1    1    0    1    1    0    0    0    1    0    1    0    1    0    0
3_1     1    1    0    1    0    0    0    1    1    1    1    1    1    0    0
3_2     1    1    0    1    0    0    0    1    1    1    1    1    1    0    0
4_1     1    0    1    0    0    1    0    1    1    0    1    1    1    1    0
4_2     1    0    1    0    0    1    0    1    1    0    1    1    1    1    0
5_1     0    1    1    1    0    0    1    1    1    0    1    1    1    1    0
5_2     0    1    1    1    0    0    1    1    1    0    1    1    1    1    0
6_1     0    0    1    0    0    0    1    1    0    0    1    0    1    0    0
6_2     0    0    1    0    0    0    1    1    0    0    1    0    1    0    0
7_1     1    1    0    0    0    0    1    0    0    0    1    1    0    0    1
7_2     1    1    0    0    0    0    1    0    0    0    1    1    0    0    1
8_1     1    0    0    0    0    0    1    1    1    1    0    0    1    1    0
8_2     1    0    0    0    0    0    1    1    1    1    0    0    1    1    0
9_1     0    0    1    1    0    1    0    1    1    0    1    0    0    0    1
9_2     0    0    1    1    0    1    0    1    1    0    1    0    0    0    1
10_1    1    1    1    0    0    0    0    1    0    1    1    1    1    0    0
10_2    1    1    1    0    0    0    0    1    0    1    1    1    1    0    0
4_65 4_66 4_67 4_68 4_69 4_70 4_71 4_72 4_73 4_74 4_75 4_76 4_77 4_78 4_79
1_1     1    0    0    1    1    0    0    1    0    1    0    0    1    0    1
1_2     1    0    0    1    1    0    0    1    0    1    0    0    1    0    1
2_1     1    1    0    0    0    1    0    0    0    0    1    0    0    0    0
2_2     1    1    0    0    0    1    0    0    0    0    1    0    0    0    0
3_1     1    0    0    0    0    1    0    0    1    1    0    0    0    0    0
3_2     1    0    0    0    0    1    0    0    1    1    0    0    0    0    0
4_1     0    1    0    1    1    1    1    1    0    1    0    0    1    1    0
4_2     0    1    0    1    1    1    1    1    0    1    0    0    1    1    0
5_1     0    0    0    0    0    0    0    1    0    0    1    0    1    0    1
5_2     0    0    0    0    0    0    0    1    0    0    1    0    1    0    1
6_1     1    1    0    1    1    0    0    1    1    0    1    1    0    0    1
6_2     1    1    0    1    1    0    0    1    1    0    1    1    0    0    1
7_1     1    1    1    0    0    1    0    0    0    0    0    0    0    0    0
7_2     1    1    1    0    0    1    0    0    0    0    0    0    0    0    0
8_1     1    1    0    0    0    1    0    0    1    0    0    1    1    1    0
8_2     1    1    0    0    0    1    0    0    1    0    0    1    1    1    0
9_1     1    1    1    0    0    0    0    0    1    1    1    0    0    0    1
9_2     1    1    1    0    0    0    0    0    1    1    1    0    0    0    1
10_1    0    1    0    0    1    1    0    0    1    1    1    1    0    0    0
10_2    0    1    0    0    1    1    0    0    1    1    1    1    0    0    0
4_80 4_81 4_82 4_83 4_84 4_85 4_86 4_87 4_88 4_89 4_90 4_91 4_92 4_93 4_94
1_1     0    0    0    1    1    1    0    1    1    0    0    1    0    1    1
1_2     0    0    0    1    1    1    0    1    1    0    0    1    0    1    1
2_1     0    0    1    1    0    0    0    1    0    1    0    1    0    0    1
2_2     0    0    1    1    0    0    0    1    0    1    0    1    0    0    1
3_1     1    0    0    1    1    0    0    0    1    1    0    1    1    0    0
3_2     1    0    0    1    1    0    0    0    1    1    0    1    1    0    0
4_1     1    0    1    1    1    0    0    0    0    0    0    1    0    0    1
4_2     1    0    1    1    1    0    0    0    0    0    0    1    0    0    1
5_1     0    1    1    0    0    1    1    1    1    1    0    1    0    1    1
5_2     0    1    1    0    0    1    1    1    1    1    0    1    0    1    1
6_1     0    1    0    0    0    1    1    0    0    0    0    1    0    0    0
6_2     0    1    0    0    0    1    1    0    0    0    0    1    0    0    0
7_1     1    0    0    0    0    1    0    1    0    0    0    1    1    1    1
7_2     1    0    0    0    0    1    0    1    0    0    0    1    1    1    1
8_1     1    1    0    0    1    1    0    0    0    0    1    1    0    0    0
8_2     1    1    0    0    1    1    0    0    0    0    1    1    0    0    0
9_1     0    0    0    0    0    1    0    0    0    0    1    1    1    1    1
9_2     0    0    0    0    0    1    0    0    0    0    1    1    1    1    1
10_1    1    1    0    1    1    1    0    1    0    1    1    0    1    0    0
10_2    1    1    0    1    1    1    0    1    0    1    1    0    1    0    0
4_95 4_96 4_97 4_98 4_99 4_100 5_1 5_2 5_3 5_4 5_5 5_6 5_7 5_8 5_9 5_10
1_1     1    1    0    1    1     1   0   0   0   1   1   1   1   0   0    1
1_2     1    1    0    1    1     1   0   0   0   1   1   1   1   0   0    1
2_1     0    0    0    1    0     1   1   0   0   0   1   0   0   1   1    0
2_2     0    0    0    1    0     1   1   0   0   0   1   0   0   1   1    0
3_1     0    1    1    1    1     0   0   0   1   0   1   1   1   0   1    1
3_2     0    1    1    1    1     0   0   0   1   0   1   1   1   0   1    1
4_1     0    0    0    0    0     0   1   0   0   1   1   1   0   0   1    0
4_2     0    0    0    0    0     0   1   0   0   1   1   1   0   0   1    0
5_1     1    0    1    0    1     0   0   0   0   0   1   0   1   1   1    1
5_2     1    0    1    0    1     0   0   0   0   0   1   0   1   1   1    1
6_1     0    1    1    0    0     1   1   1   0   0   1   0   0   1   0    1
6_2     0    1    1    0    0     1   1   1   0   0   1   0   0   1   0    1
7_1     0    0    0    0    0     1   1   1   0   1   0   1   1   1   1    0
7_2     0    0    0    0    0     1   1   1   0   1   0   1   1   1   1    0
8_1     0    0    0    0    1     1   1   1   1   0   1   0   0   1   1    0
8_2     0    0    0    0    1     1   1   1   1   0   1   0   0   1   1    0
9_1     1    1    1    0    1     1   1   1   1   1   1   1   1   1   1    0
9_2     1    1    1    0    1     1   1   1   1   1   1   1   1   1   1    0
10_1    1    0    0    1    0     0   0   0   1   1   0   1   1   1   0    1
10_2    1    0    0    1    0     0   0   0   1   1   0   1   1   1   0    1
5_11 5_12 5_13 5_14 5_15 5_16 5_17 5_18 5_19 5_20 5_21 5_22 5_23 5_24 5_25
1_1     0    1    1    1    1    0    0    1    0    1    1    1    1    0    0
1_2     0    1    1    1    1    0    0    1    0    1    1    1    1    0    0
2_1     0    1    1    1    1    0    0    0    0    0    1    0    1    1    1
2_2     0    1    1    1    1    0    0    0    0    0    1    0    1    1    1
3_1     1    0    0    1    1    0    1    1    1    0    1    1    0    1    0
3_2     1    0    0    1    1    0    1    1    1    0    1    1    0    1    0
4_1     0    1    0    0    0    1    0    1    0    0    1    1    1    1    1
4_2     0    1    0    0    0    1    0    1    0    0    1    1    1    1    1
5_1     0    0    0    1    1    0    0    0    1    0    1    1    0    1    0
5_2     0    0    0    1    1    0    0    0    1    0    1    1    0    1    0
6_1     1    1    0    0    0    1    0    1    0    1    0    1    1    1    1
6_2     1    1    0    0    0    1    0    1    0    1    0    1    1    1    1
7_1     1    0    1    0    0    0    1    0    1    1    0    1    0    1    1
7_2     1    0    1    0    0    0    1    0    1    1    0    1    0    1    1
8_1     0    1    0    1    1    0    1    0    0    1    0    0    0    1    0
8_2     0    1    0    1    1    0    1    0    0    1    0    0    0    1    0
9_1     1    0    1    0    1    0    1    0    1    0    0    0    0    1    0
9_2     1    0    1    0    1    0    1    0    1    0    0    0    0    1    0
10_1    0    0    1    1    1    0    0    0    0    1    0    1    1    0    0
10_2    0    0    1    1    1    0    0    0    0    1    0    1    1    0    0
5_26 5_27 5_28 5_29 5_30 5_31 5_32 5_33 5_34 5_35 5_36 5_37 5_38 5_39 5_40
1_1     0    1    0    0    0    0    1    0    1    1    0    1    0    0    1
1_2     0    1    0    0    0    0    1    0    1    1    0    1    0    0    1
2_1     0    1    1    0    0    1    1    0    1    1    0    0    1    0    0
2_2     0    1    1    0    0    1    1    0    1    1    0    0    1    0    0
3_1     0    1    0    0    1    1    1    1    0    0    0    1    1    0    1
3_2     0    1    0    0    1    1    1    1    0    0    0    1    1    0    1
4_1     1    1    1    1    1    1    1    0    0    0    0    1    1    1    0
4_2     1    1    1    1    1    1    1    0    0    0    0    1    1    1    0
5_1     1    1    0    0    1    0    0    0    0    1    1    0    1    1    1
5_2     1    1    0    0    1    0    0    0    0    1    1    0    1    1    1
6_1     1    1    1    0    1    1    0    0    0    0    0    1    1    0    1
6_2     1    1    1    0    1    1    0    0    0    0    0    1    1    0    1
7_1     0    0    1    1    1    0    0    1    0    0    0    1    0    1    0
7_2     0    0    1    1    1    0    0    1    0    0    0    1    0    1    0
8_1     0    1    0    0    0    0    0    1    1    0    1    0    0    0    0
8_2     0    1    0    0    0    0    0    1    1    0    1    0    0    0    0
9_1     1    1    0    1    1    1    0    1    1    1    1    0    1    1    1
9_2     1    1    0    1    1    1    0    1    1    1    1    0    1    1    1
10_1    1    1    1    0    1    0    0    1    0    1    1    0    0    1    1
10_2    1    1    1    0    1    0    0    1    0    1    1    0    0    1    1
5_41 5_42 5_43 5_44 5_45 5_46 5_47 5_48 5_49 5_50 5_51 5_52 5_53 5_54 5_55
1_1     0    0    0    1    0    1    1    0    1    0    1    0    0    0    0
1_2     0    0    0    1    0    1    1    0    1    0    1    0    0    0    0
2_1     0    0    0    0    1    1    1    1    0    0    0    1    0    1    1
2_2     0    0    0    0    1    1    1    1    0    0    0    1    0    1    1
3_1     0    0    0    0    0    1    0    0    0    1    1    0    0    1    1
3_2     0    0    0    0    0    1    0    0    0    1    1    0    0    1    1
4_1     1    1    0    0    1    0    0    0    0    0    1    1    0    1    0
4_2     1    1    0    0    1    0    0    0    0    0    1    1    0    1    0
5_1     1    1    0    1    1    1    1    1    0    0    1    1    1    1    0
5_2     1    1    0    1    1    1    1    1    0    0    1    1    1    1    0
6_1     1    0    1    0    1    0    1    0    0    1    0    0    0    0    0
6_2     1    0    1    0    1    0    1    0    0    1    0    0    0    0    0
7_1     0    1    0    0    0    1    1    0    1    1    0    1    0    1    1
7_2     0    1    0    0    0    1    1    0    1    1    0    1    0    1    1
8_1     1    0    0    1    1    1    0    0    1    1    1    0    1    0    1
8_2     1    0    0    1    1    1    0    0    1    1    1    0    1    0    1
9_1     1    0    1    0    1    0    0    0    1    0    1    0    0    0    0
9_2     1    0    1    0    1    0    0    0    1    0    1    0    0    0    0
10_1    1    0    1    0    1    0    1    0    0    1    1    1    0    1    1
10_2    1    0    1    0    1    0    1    0    0    1    1    1    0    1    1
5_56 5_57 5_58 5_59 5_60 5_61 5_62 5_63 5_64 5_65 5_66 5_67 5_68 5_69 5_70
1_1     0    1    1    1    1    1    1    1    0    1    0    0    0    0    0
1_2     0    1    1    1    1    1    1    1    0    1    0    0    0    0    0
2_1     1    0    1    0    1    1    0    0    0    1    1    0    1    0    1
2_2     1    0    1    0    1    1    0    0    0    1    1    0    1    0    1
3_1     1    1    0    0    0    1    1    0    0    0    1    1    0    1    1
3_2     1    1    0    0    0    1    1    0    0    0    1    1    0    1    1
4_1     0    1    0    0    1    0    1    0    1    1    0    1    1    1    1
4_2     0    1    0    0    1    0    1    0    1    1    0    1    1    1    1
5_1     0    0    0    1    0    0    1    1    1    0    1    1    1    0    0
5_2     0    0    0    1    0    0    1    1    1    0    1    1    1    0    0
6_1     0    1    0    1    1    0    1    1    1    0    0    0    1    1    1
6_2     0    1    0    1    1    0    1    1    1    0    0    0    1    1    1
7_1     0    1    0    1    0    0    0    0    0    1    1    1    0    1    1
7_2     0    1    0    1    0    0    0    0    0    1    1    1    0    1    1
8_1     0    0    0    1    0    0    0    0    1    1    0    0    0    0    1
8_2     0    0    0    1    0    0    0    0    1    1    0    0    0    0    1
9_1     1    1    1    0    1    0    0    1    1    1    0    0    1    0    0
9_2     1    1    1    0    1    0    0    1    1    1    0    0    1    0    0
10_1    1    1    0    0    1    0    0    0    1    0    1    1    1    0    1
10_2    1    1    0    0    1    0    0    0    1    0    1    1    1    0    1
5_71 5_72 5_73 5_74 5_75 5_76 5_77 5_78 5_79 5_80 5_81 5_82 5_83 5_84 5_85
1_1     0    1    0    1    0    1    0    0    1    0    0    0    1    0    1
1_2     0    1    0    1    0    1    0    0    1    0    0    0    1    0    1
2_1     1    0    0    0    1    1    1    1    0    1    0    0    0    1    1
2_2     1    0    0    0    1    1    1    1    0    1    0    0    0    1    1
3_1     1    1    0    1    0    1    1    0    1    0    1    0    0    1    1
3_2     1    1    0    1    0    1    1    0    1    0    1    0    0    1    1
4_1     1    0    1    0    1    1    1    1    1    0    1    1    0    0    1
4_2     1    0    1    0    1    1    1    1    1    0    1    1    0    0    1
5_1     0    0    0    0    0    0    0    0    1    0    1    1    0    1    0
5_2     0    0    0    0    0    0    0    0    1    0    1    1    0    1    0
6_1     1    0    0    0    0    1    0    0    0    0    0    0    1    1    1
6_2     1    0    0    0    0    1    0    0    0    0    0    0    1    1    1
7_1     0    1    1    0    0    0    1    1    0    1    0    1    0    1    0
7_2     0    1    1    0    0    0    1    1    0    1    0    1    0    1    0
8_1     1    1    0    1    0    1    0    0    1    0    1    1    0    1    0
8_2     1    1    0    1    0    1    0    0    1    0    1    1    0    1    0
9_1     1    1    1    1    0    1    1    1    1    0    1    1    1    1    0
9_2     1    1    1    1    0    1    1    1    1    0    1    1    1    1    0
10_1    1    1    0    1    1    1    0    1    1    0    0    1    0    0    0
10_2    1    1    0    1    1    1    0    1    1    0    0    1    0    0    0
5_86 5_87 5_88 5_89 5_90 5_91 5_92 5_93 5_94 5_95 5_96 5_97 5_98 5_99
1_1     0    1    1    1    1    0    0    0    0    1    1    0    1    0
1_2     0    1    1    1    1    0    0    0    0    1    1    0    1    0
2_1     1    1    0    1    1    0    0    0    1    0    1    0    0    0
2_2     1    1    0    1    1    0    0    0    1    0    1    0    0    0
3_1     0    1    0    1    1    1    1    0    0    1    1    0    1    0
3_2     0    1    0    1    1    1    1    0    0    1    1    0    1    0
4_1     0    0    0    1    0    0    1    1    0    0    0    1    0    0
4_2     0    0    0    1    0    0    1    1    0    0    0    1    0    0
5_1     0    0    0    1    1    0    1    0    1    0    1    0    0    0
5_2     0    0    0    1    1    0    1    0    1    0    1    0    0    0
6_1     0    0    0    0    1    0    0    0    0    1    0    0    1    1
6_2     0    0    0    0    1    0    0    0    0    1    0    0    1    1
7_1     1    1    1    0    0    1    1    0    0    1    1    1    0    1
7_2     1    1    1    0    0    1    1    0    0    1    1    1    0    1
8_1     1    0    0    1    1    0    0    0    0    0    0    0    1    1
8_2     1    0    0    1    1    0    0    0    0    0    0    0    1    1
9_1     1    1    1    1    1    0    0    1    1    1    1    1    0    1
9_2     1    1    1    1    1    0    0    1    1    1    1    1    0    1
10_1    0    1    1    0    0    1    0    1    1    1    0    0    1    1
10_2    0    1    1    0    0    1    0    1    1    1    0    0    1    1
5_100
1_1      1
1_2      1
2_1      0
2_2      0
3_1      0
3_2      0
4_1      0
4_2      0
5_1      0
5_2      0
6_1      0
6_2      0
7_1      0
7_2      0
8_1      1
8_2      1
9_1      1
9_2      1
10_1     0
10_2     0
> dim(basePopHaplo)
[1]  20 500
> class(basePopHaplo)
[1] "matrix" "array"
> str(basePopHaplo)
int [1:20, 1:500] 0 0 1 1 0 0 1 1 0 0 ...
- attr(*, "dimnames")=List of 2
..$ : chr [1:20] "1_1" "1_2" "2_1" "2_2" ...
..$ : chr [1:500] "1_1" "1_2" "1_3" "1_4" ...
> founderGenome=quickHaplo(nInd =5, nChr = 5, segSites = 10, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=10, mean=5, var=1)
basePop=newPop(founderGenome)
> founderGenome=quickHaplo(nInd = 10, nChr = 1, segSites = 10, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=10, mean=10, var=1)
basePop=newPop(founderGenome)
> str(basePop)
Formal class 'Pop' [package "AlphaSimR"] with 18 slots
..@ id     : chr [1:10] "1" "2" "3" "4" ...
..@ iid    : int [1:10] 1 2 3 4 5 6 7 8 9 10
..@ mother : chr [1:10] "0" "0" "0" "0" ...
..@ father : chr [1:10] "0" "0" "0" "0" ...
..@ sex    : chr [1:10] "H" "H" "H" "H" ...
..@ nTraits: int 1
..@ gv     : num [1:10, 1] 11.69 10.21 9.89 10.3 10.41 ...
.. ..- attr(*, "dimnames")=List of 2
.. .. ..$ : NULL
.. .. ..$ : chr "Trait1"
..@ pheno  : num [1:10, 1] NA NA NA NA NA NA NA NA NA NA
.. ..- attr(*, "dimnames")=List of 2
.. .. ..$ : NULL
.. .. ..$ : chr "Trait1"
..@ ebv    : num[1:10, 0 ]
..@ gxe    :List of 1
.. ..$ : NULL
..@ fixEff : int [1:10] 1 1 1 1 1 1 1 1 1 1
..@ misc   :List of 10
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
..@ miscPop: list()
..@ nInd   : int 10
..@ nChr   : int 1
..@ ploidy : int 2
..@ nLoci  : int 10
..@ geno   :List of 1
.. ..$ : raw [1:2, 1:2, 1:10] a5 da a5 da ...
> dim(basePop)
NULL
> class(basePop)
[1] "Pop"
attr(,"package")
[1] "AlphaSimR"
> basePopHaplo=pullSegSiteHaplo(basePop)
> dim(basePopHaplo)
[1] 20 10
> class(basePopHaplo)
[1] "matrix" "array"
> str(basePopHaplo)
int [1:20, 1:10] 1 1 0 0 1 1 0 0 1 1 ...
- attr(*, "dimnames")=List of 2
..$ : chr [1:20] "1_1" "1_2" "2_1" "2_2" ...
..$ : chr [1:10] "1_1" "1_2" "1_3" "1_4" ...
> founderGenome=quickHaplo(nInd = 10, nChr = 1, segSites = 5, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=10, mean=10, var=1)
basePop=newPop(founderGenome)
Error in private$.pickLoci(nQtlPerChr) :
  sapply(pot, length) >= nSitesPerChr is not TRUE
> founderGenome=quickHaplo(nInd = 10, nChr = 1, segSites = 5, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=5, mean=10, var=1)
basePop=newPop(founderGenome)
> basePopHaplo=pullSegSiteHaplo(basePop)
> dim(basePopHaplo)
[1] 20  5
> founderGenome=quickHaplo(nInd = 9, nChr = 1, segSites = 5, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=5, mean=10, var=1)
basePop=newPop(founderGenome)
basePopHaplo=pullSegSiteHaplo(basePop)
dim(basePopHaplo)
[1] 18  5
> basePopHaplo$
  basePopHaplo$
  > basePopHaplo@
  basePopHaplo@
  > basePopHaplo@
  basePopHaplo@
  > basePopHaplo
1_1 1_2 1_3 1_4 1_5
1_1   1   0   0   0   1
1_2   1   0   0   0   1
2_1   0   0   0   1   1
2_2   0   0   0   1   1
3_1   1   0   0   1   1
3_2   1   0   0   1   1
4_1   0   0   1   0   1
4_2   0   0   1   0   1
5_1   0   1   0   0   1
5_2   0   1   0   0   1
6_1   0   1   0   1   1
6_2   0   1   0   1   1
7_1   0   1   0   0   1
7_2   0   1   0   0   1
8_1   1   1   1   0   1
8_2   1   1   1   0   1
9_1   0   0   1   0   0
9_2   0   0   1   0   0
> ls()
[1] "basePop"       "basePopHaplo"  "data_path"     "founderGenome"
[5] "geno"          "SP"
> dim(geno)
[1]    344 207981
> geno[1:5,1:5]
m1 m2 m3 m4 m5
CS76161_1  0  0  0  0  0
CS76189_1  0  0  0  0  0
CS76160_1  0  2  0  0  0
CS76253_1  0  0  2  0  0
CS76158_1  0  0  0  0  0
> pullSnpGeno(basePop)
Error in simParam$snpChips[[snpChip]] : subscript out of bounds
> pullSnpGeno(basePop)
basePopGeno=pullSegSiteGeno(basePop)
> basePopGeno [, 1:5]
1_1 1_2 1_3 1_4 1_5
1   2   0   0   0   2
2   0   0   0   2   2
3   2   0   0   2   2
4   0   0   2   0   2
5   0   2   0   0   2
6   0   2   0   2   2
7   0   2   0   0   2
8   2   2   2   0   2
9   0   0   2   0   0
> dim(basePopGeno)
[1] 9 5
> class(basePopGeno)
[1] "matrix" "array"
> str(basePopGeno)
int [1:9, 1:5] 2 0 2 0 0 0 0 2 0 0 ...
- attr(*, "dimnames")=List of 2
..$ : chr [1:9] "1" "2" "3" "4" ...
..$ : chr [1:5] "1_1" "1_2" "1_3" "1_4" ...
> rowSums(basePopHaplo)
1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2 9_1 9_2
2   2   2   2   3   3   2   2   2   2   3   3   2   2   4   4   1   1

cross12=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "11" "1"    "2"
[2,] "12" "1"    "2"
[3,] "13" "1"    "2"
[4,] "14" "1"    "2"
[5,] "15" "1"    "2"
[6,] "16" "1"    "2"
cross13=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "11" "1"    "2"
[2,] "12" "1"    "2"
[3,] "13" "1"    "2"
[4,] "14" "1"    "2"
[5,] "15" "1"    "2"
[6,] "16" "1"    "2"
cross23=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "31" "1"    "2"
[2,] "32" "1"    "2"
[3,] "33" "1"    "2"
[4,] "34" "1"    "2"
[5,] "35" "1"    "2"
[6,] "36" "1"    "2"
hist(x(cross12, cross13, cross23), 
     freq = NULL,
     xlim = range(breaks=seq(from=range(1), to=range(2),)
     ylab="Number of mutations",
     abline(v=nPar1, col="blue", lwd=3),
     abline(v=mPar2, col="red", lwd=3),
     abline(v=nPar3, col="green",lwd=3),
     plot = TRUE)
######################################################################################################################################################

dim(geno)
[1]    344 207981

geno[1:5,1:5]
m1 m2 m3 m4 m5
CS76161_1  0  0  0  0  0
CS76189_1  0  0  0  0  0
CS76160_1  0  2  0  0  0
CS76253_1  0  0  2  0  0
CS76158_1  0  0  0  0  0

#Create Haplotypes from Founder Population
founderGenome=quickHaplo(nInd = 10, nChr = 1, segSites = 5, inbred=TRUE)
range
#Set Parameters
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=5, mean=10, var=1)

#Create Base Population
basePop=newPop(founderGenome)
basePopHaplo=pullSegSiteHaplo(basePop)

dim(basePopHaplo)
basePopHaplo=pullSegSiteHaplo(basePop)
basePopHaplo [, 1:5]
1_1 1_2 1_3 1_4 1_5
1_1    1   1   1   1   1
1_2    1   1   1   1   1
2_1    0   0   0   1   1
2_2    0   0   0   1   1
3_1    0   0   1   1   0
3_2    0   0   1   1   0
4_1    0   0   0   0   0
4_2    0   0   0   0   0
5_1    1   1   0   0   1
5_2    1   1   0   0   1
6_1    0   1   1   0   1
6_2    0   1   1   0   1
7_1    1   1   1   1   1
7_2    1   1   1   1   1
8_1    1   1   1   1   1
8_2    1   1   1   1   1
9_1    0   1   1   0   1
9_2    0   1   1   0   1
10_1   0   1   0   1   0
10_2   0   1   0   1   0

basePopHaplo

basePopGeno=pullSegSiteGeno(basePop)
basePopGeno [, 1:5]
rowSums(basePopHaplo)
rowSums(basePopGeno)
////////////////////////////////////////////////////////////////////////////////////////////////
1.

cross12=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "10" "1"    "2"
[2,] "11" "1"    "2"
[3,] "12" "1"    "2"
[4,] "13" "1"    "2"
[5,] "14" "1"    "2"
[6,] "15" "1"    "2"

cross12Geno=pullSegSiteGeno(cross12)
pullSegSiteGeno(cross12)
1_1 1_2 1_3 1_4 1_5
10   1   2   1   2   1
11   1   2   1   2   1
12   1   2   1   2   1
13   1   2   1   2   1
14   1   2   1   2   1
15   1   2   1   2   1
16   1   2   1   2   1
17   1   2   1   2   1
18   1   2   1   2   1
19   1   2   1   2   1

ncross12=rowSums(cross12Geno)
rowSums(cross12Geno)
10 11 12 13 14 15 16 17 18 19
7  7  7  7  7  7  7  7  7  7

cross13=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "10" "1"    "2"
[2,] "11" "1"    "2"
[3,] "12" "1"    "2"
[4,] "13" "1"    "2"
[5,] "14" "1"    "2"
[6,] "15" "1"    "2"

cross13Geno=pullSegSiteGeno(cross13)
pullSegSiteGeno(cross13)
1_1 1_2 1_3 1_4 1_5
20   1   2   1   2   1
21   1   2   1   2   1
22   1   2   1   2   1
23   1   2   1   2   1
24   1   2   1   2   1
25   1   2   1   2   1
26   1   2   1   2   1
27   1   2   1   2   1
28   1   2   1   2   1
29   1   2   1   2   1

ncross13=rowSums(cross13Geno)
rowSums(cross13Geno)
20 21 22 23 24 25 26 27 28 29
7  7  7  7  7  7  7  7  7  7
/////////////////////////////////////////////////////////////////////////////////////////

cross23=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "10" "1"    "2"
[2,] "11" "1"    "2"
[3,] "12" "1"    "2"
[4,] "13" "1"    "2"
[5,] "14" "1"    "2"
[6,] "15" "1"    "2"

cross23Geno=pullSegSiteGeno(cross23)
pullSegSiteGeno(cross23)
1_1 1_2 1_3 1_4 1_5
30   1   2   1   2   1
31   1   2   1   2   1
32   1   2   1   2   1
33   1   2   1   2   1
34   1   2   1   2   1
35   1   2   1   2   1
36   1   2   1   2   1
37   1   2   1   2   1
38   1   2   1   2   1
39   1   2   1   2   1

ncross23=rowSums(cross23Geno)
rowSums(cross23Geno)
30 31 32 33 34 35 36 37 38 39
7  7  7  7  7  7  7  7  7  7

#######################################################################################################################################################
dim(geno)
[1]    344 207981
> geno[1:5,1:5]
m1 m2 m3 m4 m5
CS76161_1  0  0  0  0  0
CS76189_1  0  0  0  0  0
CS76160_1  0  2  0  0  0
CS76253_1  0  0  2  0  0
CS76158_1  0  0  0  0  0
> founderGenome=quickHaplo(nInd = 9, nChr = 1, segSites = 5, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=5, mean=10, var=1)
basePop=newPop(founderGenome)
> basePopHaplo=pullSegSiteHaplo(basePop)
Error: unexpected '>' in ">"
> founderGenome=quickHaplo(nInd = 9, nChr = 1, segSites = 5, inbred=TRUE)
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=5, mean=10, var=1)
> basePop=newPop(founderGenome)
basePopHaplo=pullSegSiteHaplo(basePop)
dim(basePopHaplo)
[1] 18  5
> basePopHaplo=pullSegSiteHaplo(basePop)
basePopHaplo [, 1:5]
basePopHaplo
1_1 1_2 1_3 1_4 1_5
1_1   1   1   1   1   1
1_2   1   1   1   1   1
2_1   0   1   0   1   0
2_2   0   1   0   1   0
3_1   0   0   1   0   0
3_2   0   0   1   0   0
4_1   0   0   0   1   0
4_2   0   0   0   1   0
5_1   1   0   0   1   1
5_2   1   0   0   1   1
6_1   0   0   1   1   1
6_2   0   0   1   1   1
7_1   1   1   1   1   1
7_2   1   1   1   1   1
8_1   0   0   0   1   0
8_2   0   0   0   1   0
9_1   1   0   1   0   1
9_2   1   0   1   0   1
1_1 1_2 1_3 1_4 1_5
1_1   1   1   1   1   1
1_2   1   1   1   1   1
2_1   0   1   0   1   0
2_2   0   1   0   1   0
3_1   0   0   1   0   0
3_2   0   0   1   0   0
4_1   0   0   0   1   0
4_2   0   0   0   1   0
5_1   1   0   0   1   1
5_2   1   0   0   1   1
6_1   0   0   1   1   1
6_2   0   0   1   1   1
7_1   1   1   1   1   1
7_2   1   1   1   1   1
8_1   0   0   0   1   0
8_2   0   0   0   1   0
9_1   1   0   1   0   1
9_2   1   0   1   0   1
> basePopGeno=pullSegSiteGeno(basePop)
basePopGeno [, 1:5]
1_1 1_2 1_3 1_4 1_5
1   2   2   2   2   2
2   0   2   0   2   0
3   0   0   2   0   0
4   0   0   0   2   0
5   2   0   0   2   2
6   0   0   2   2   2
7   2   2   2   2   2
8   0   0   0   2   0
9   2   0   2   0   2
> rowSums(basePopHaplo)
1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2 9_1 9_2
5   5   2   2   1   1   1   1   3   3   3   3   5   5   1   1   3   3
> rowSums(basePopGeno)
1  2  3  4  5  6  7  8  9
10  4  2  2  6  6 10  2  6
> cross12=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
> head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "10" "1"    "2"
[2,] "11" "1"    "2"
[3,] "12" "1"    "2"
[4,] "13" "1"    "2"
[5,] "14" "1"    "2"
[6,] "15" "1"    "2"
> cross12Geno=pullSegSiteGeno(cross12)
> pullSegSiteGeno(cross12)
1_1 1_2 1_3 1_4 1_5
10   1   2   1   2   1
11   1   2   1   2   1
12   1   2   1   2   1
13   1   2   1   2   1
14   1   2   1   2   1
15   1   2   1   2   1
16   1   2   1   2   1
17   1   2   1   2   1
18   1   2   1   2   1
19   1   2   1   2   1
> ncross12=rowSums(cross12Geno)
> rowSums(cross12Geno)
10 11 12 13 14 15 16 17 18 19
7  7  7  7  7  7  7  7  7  7
> cross13=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
> head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "10" "1"    "2"
[2,] "11" "1"    "2"
[3,] "12" "1"    "2"
[4,] "13" "1"    "2"
[5,] "14" "1"    "2"
[6,] "15" "1"    "2"
> cross13Geno=pullSegSiteGeno(cross13)
> pullSegSiteGeno(cross13)
1_1 1_2 1_3 1_4 1_5
20   1   2   1   2   1
21   1   2   1   2   1
22   1   2   1   2   1
23   1   2   1   2   1
24   1   2   1   2   1
25   1   2   1   2   1
26   1   2   1   2   1
27   1   2   1   2   1
28   1   2   1   2   1
29   1   2   1   2   1
> ncross13=rowSums(cross13Geno)
> rowSums(cross13Geno)
20 21 22 23 24 25 26 27 28 29
7  7  7  7  7  7  7  7  7  7
> cross23=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
> head(cbind(id=cross12@id, mother=cross12@mother, father=cross12@father))
id   mother father
[1,] "10" "1"    "2"
[2,] "11" "1"    "2"
[3,] "12" "1"    "2"
[4,] "13" "1"    "2"
[5,] "14" "1"    "2"
[6,] "15" "1"    "2"
> cross23Geno=pullSegSiteGeno(cross23)
> pullSegSiteGeno(cross23)
1_1 1_2 1_3 1_4 1_5
30   1   2   1   2   1
31   1   2   1   2   1
32   1   2   1   2   1
33   1   2   1   2   1
34   1   2   1   2   1
35   1   2   1   2   1
36   1   2   1   2   1
37   1   2   1   2   1
38   1   2   1   2   1
39   1   2   1   2   1
> ncross23=rowSums(cross23Geno)
> rowSums(cross23Geno)
30 31 32 33 34 35 36 37 38 39
7  7  7  7  7  7  7  7  7  7
>
  > #Create Haplotypes from Founder Population
  founderGenome=quickHaplo(nInd = 10, nChr = 1, segSites = 5, inbred=TRUE)

#Set Parameters
SP=SimParam$new(founderGenome)
SP$addTraitA(nQtlPerChr=5, mean=10, var=1)

#Create Base Population
basePop=newPop(founderGenome)
basePopHaplo=pullSegSiteHaplo(basePop)

dim(basePopHaplo)
basePopHaplo=pullSegSiteHaplo(basePop)
basePopHaplo [, 1:5]
basePopHaplo

basePopGeno=pullSegSiteGeno(basePop)
basePopGeno [, 1:5]
rowSums(basePopHaplo)
rowSums(basePopGeno)


cross12=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 10)
[1] 20  5
1_1 1_2 1_3 1_4 1_5
1_1    0   1   1   1   0
1_2    0   1   1   1   0
2_1    0   0   0   0   0
2_2    0   0   0   0   0
3_1    1   1   0   0   0
3_2    1   1   0   0   0
4_1    0   1   0   1   0
4_2    0   1   0   1   0
5_1    0   1   0   0   0
5_2    0   1   0   0   0
6_1    0   0   1   1   1
6_2    0   0   1   1   1
7_1    1   0   0   0   0
7_2    1   0   0   0   0
8_1    1   0   1   0   0
8_2    1   0   1   0   0
9_1    1   0   1   1   0
9_2    1   0   1   1   0
10_1   1   0   0   1   0
10_2   1   0   0   1   0
1_1 1_2 1_3 1_4 1_5
1_1    0   1   1   1   0
1_2    0   1   1   1   0
2_1    0   0   0   0   0
2_2    0   0   0   0   0
3_1    1   1   0   0   0
3_2    1   1   0   0   0
4_1    0   1   0   1   0
4_2    0   1   0   1   0
5_1    0   1   0   0   0
5_2    0   1   0   0   0
6_1    0   0   1   1   1
6_2    0   0   1   1   1
7_1    1   0   0   0   0
7_2    1   0   0   0   0
8_1    1   0   1   0   0
8_2    1   0   1   0   0
9_1    1   0   1   1   0
9_2    1   0   1   1   0
10_1   1   0   0   1   0
10_2   1   0   0   1   0
1_1 1_2 1_3 1_4 1_5
1    0   2   2   2   0
2    0   0   0   0   0
3    2   2   0   0   0
4    0   2   0   2   0
5    0   2   0   0   0
6    0   0   2   2   2
7    2   0   0   0   0
8    2   0   2   0   0
9    2   0   2   2   0
10   2   0   0   2   0
1_1  1_2  2_1  2_2  3_1  3_2  4_1  4_2  5_1  5_2  6_1  6_2  7_1  7_2  8_1  8_2
3    3    0    0    2    2    2    2    1    1    3    3    1    1    2    2
9_1  9_2 10_1 10_2
3    3    2    2
1  2  3  4  5  6  7  8  9 10
6  0  4  4  2  6  2  4  6  4
>
  >
  >
  >
  > cross12
An object of class "Pop"
Ploidy: 2
Individuals: 10
Chromosomes: 1
Loci: 5
Traits: 1
> cross12=makeCross(pop = basePop, crossPlan = matrix(c(1,2), ncol = 2), nProgeny = 15)
> cross12
An object of class "Pop"
Ploidy: 2
Individuals: 15
Chromosomes: 1
Loci: 5
Traits: 1
> basePopGeno
1_1 1_2 1_3 1_4 1_5
1    0   2   2   2   0
2    0   0   0   0   0
3    2   2   0   0   0
4    0   2   0   2   0
5    0   2   0   0   0
6    0   0   2   2   2
7    2   0   0   0   0
8    2   0   2   0   0
9    2   0   2   2   0
10   2   0   0   2   0
> dim(cross12)
NULL
> class(cross12)
[1] "Pop"
attr(,"package")
[1] "AlphaSimR"
> str(cross12)
Formal class 'Pop' [package "AlphaSimR"] with 18 slots
..@ id     : chr [1:15] "21" "22" "23" "24" ...
..@ iid    : int [1:15] 21 22 23 24 25 26 27 28 29 30 ...
..@ mother : chr [1:15] "1" "1" "1" "1" ...
..@ father : chr [1:15] "2" "2" "2" "2" ...
..@ sex    : chr [1:15] "H" "H" "H" "H" ...
..@ nTraits: int 1
..@ gv     : num [1:15, 1] 10.4 10.4 10.4 10.4 10.4 ...
.. ..- attr(*, "dimnames")=List of 2
.. .. ..$ : NULL
.. .. ..$ : chr "Trait1"
..@ pheno  : num [1:15, 1] NA NA NA NA NA NA NA NA NA NA ...
.. ..- attr(*, "dimnames")=List of 2
.. .. ..$ : NULL
.. .. ..$ : chr "Trait1"
..@ ebv    : num[1:15, 0 ]
..@ gxe    :List of 1
.. ..$ : NULL
..@ fixEff : int [1:15] 1 1 1 1 1 1 1 1 1 1 ...
..@ misc   :List of 15
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
.. ..$ : NULL
..@ miscPop: list()
..@ nInd   : int 15
..@ nChr   : int 1
..@ ploidy : int 2
..@ nLoci  : int 5
..@ geno   :List of 1
.. ..$ : raw [1, 1:2, 1:15] ee 60 ee 60 ...
> cross12$
  cross12$
  > cross12@
  cross12@id       cross12@father   cross12@gv       cross12@gxe      cross12@miscPop  cross12@ploidy
cross12@iid      cross12@sex      cross12@pheno    cross12@fixEff   cross12@nInd     cross12@nLoci
cross12@mother   cross12@nTraits  cross12@ebv      cross12@misc     cross12@nChr     cross12@geno
> cross12@id
[1] "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31" "32" "33" "34" "35"
> cross12@iid
[1] 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
> cross12@mother
[1] "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1" "1"
> cross12@father
[1] "2" "2" "2" "2" "2" "2" "2" "2" "2" "2" "2" "2" "2" "2" "2"
> cross12@sex
[1] "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H" "H"
> cross12@nTraits
[1] 1
> cross12@gv
Trait1
[1,] 10.40813
[2,] 10.40813
[3,] 10.40813
[4,] 10.40813
[5,] 10.40813
[6,] 10.40813
[7,] 10.40813
[8,] 10.40813
[9,] 10.40813
[10,] 10.40813
[11,] 10.40813
[12,] 10.40813
[13,] 10.40813
[14,] 10.40813
[15,] 10.40813
> basePop
basePop       basePopGeno   basePopHaplo
> basePopGeno
1_1 1_2 1_3 1_4 1_5
1    0   2   2   2   0
2    0   0   0   0   0
3    2   2   0   0   0
4    0   2   0   2   0
5    0   2   0   0   0
6    0   0   2   2   2
7    2   0   0   0   0
8    2   0   2   0   0
9    2   0   2   2   0
10   2   0   0   2   0
> cross12@pheno
Trait1
[1,]     NA
[2,]     NA
[3,]     NA
[4,]     NA
[5,]     NA
[6,]     NA
[7,]     NA
[8,]     NA
[9,]     NA
[10,]     NA
[11,]     NA
[12,]     NA
[13,]     NA
[14,]     NA
[15,]     NA
> cross12@ebv

[1,]
[2,]
[3,]
[4,]
[5,]
[6,]
[7,]
[8,]
[9,]
[10,]
[11,]
[12,]
[13,]
[14,]
[15,]
> cross12@geno
[[1]]
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60


> class(cross12@geno)
[1] "list"
> length(cross12@geno)
[1] 1
> length(cross12@geno[1])
[1] 1
> length(cross12@geno[[1]])
[1] 30
> length(cross12@geno[[2]])
Error in cross12@geno[[2]] : subscript out of bounds
> basePopGeno
1_1 1_2 1_3 1_4 1_5
1    0   2   2   2   0
2    0   0   0   0   0
3    2   2   0   0   0
4    0   2   0   2   0
5    0   2   0   0   0
6    0   0   2   2   2
7    2   0   0   0   0
8    2   0   2   0   0
9    2   0   2   2   0
10   2   0   0   2   0
> unlist(cross12@geno[[1]])
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60

> unlist(unlist(cross12@geno[[1]]))
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60

> unlist(unlist(unlist(cross12@geno[[1]])))
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60

> cross12@geno[[1]][1]
[1] ee
> cross12@geno[[1]]
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60

> cross12@geno[[1]][1:222]
[1] ee 60 ee 60 ee 60 ee 60 ee 60 ee 60 ee 60 ee 60 ee 60 ee 60 ee 60 ee 60 ee
[26] 60 ee 60 ee 60 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
[51] 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
[76] 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
[101] 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
[126] 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
[151] 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
[176] 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
[201] 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00
> cross12@geno[[1]]
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60

> cross12@geno[1]

[[1]]
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60


> cross12@geno[[1]]

, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60

> cross12@geno[1]

[[1]]
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60


> cross12@geno[1][1]

[[1]]
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60


> cross12@geno[1][1][1]

[[1]]
, , 1

[,1] [,2]
[1,]   ee   60

, , 2

[,1] [,2]
[1,]   ee   60

, , 3

[,1] [,2]
[1,]   ee   60

, , 4

[,1] [,2]
[1,]   ee   60

, , 5

[,1] [,2]
[1,]   ee   60

, , 6

[,1] [,2]
[1,]   ee   60

, , 7

[,1] [,2]
[1,]   ee   60

, , 8

[,1] [,2]
[1,]   ee   60

, , 9

[,1] [,2]
[1,]   ee   60

, , 10

[,1] [,2]
[1,]   ee   60

, , 11

[,1] [,2]
[1,]   ee   60

, , 12

[,1] [,2]
[1,]   ee   60

, , 13

[,1] [,2]
[1,]   ee   60

, , 14

[,1] [,2]
[1,]   ee   60

, , 15

[,1] [,2]
[1,]   ee   60


> length(cross12@geno[1])

[1] 1

################################################################################
#combine genotypes of all individuals in one object
allGeno=rbind(basePopGeno, cross12Geno)
hist(cross12, xlin=rangeVal, breaks=bins, ylab="Genetic values")
abline(v=gvPar1, col="blue", lwd=3), 
abline(v=gvPar2, col="red", lwd=3),
abline(v=gvPar1+gvPar2/2, col="black", lwd=3, lty=2))


  

















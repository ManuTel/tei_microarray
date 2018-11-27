# Differential gene expression microarray analysis in paternal diet induced Trangenerational epigenetic inheritance model in Drosophila melanogaster
All the microarray analysis related to Transgenerational epigenetic inheritance in drosophila melanogaster.

The project deals with the gene expression analysis in Drosophila model of Transgenerational epigenetic inheritance. Tha main analysis consist of following differential gene expression profile

    - **Case**: In F1 and F2 generation, Case means flies having High sugar diet fed paternal ancestry
    - **Control**: means flies having control diet fed paternal ancestry.

| Sr.no.     | Generation   | Gender | Comparision      | Diet       |
| ---------- |:------------:|:------:|:----------------:|:---------: |
| 01         | F0           | Male   | High sugar vs CD | *NA*       |
| 02         | F0           | Female | High sugar vs CD | *NA*       |
| 03         | F1           | Male   | Case vs Control  | Control    |
| 04         | F1           | Male   | Case vs Control  | High sugar |
| 05         | F1           | Female | Case vs Control  | Control    |
| 06         | F1           | Female | Case vs Control  | High sugar |
| 07         | F2           | Male   | Case vs Control  | Control    |
| 08         | F2           | Male   | Case vs Control  | High sugar |
| 09         | F2           | Female | Case vs Control  | Control    |
| 10         | F2           | Female | Case vs Control  | High sugar |

## Objectives of the analysis is
1. To find out the effect of 2-week high sugar diet on the gene expression profile in both sexes of *Drosophila melanogaster*
2. F1 generation:
    - What is the effect of F0 male parent diet(high sugar vs control) on the F1 progeny's gene expression pattern.
    - Is it F1 progeny sex-specific?
    - what fraction of genes are expressed in same direction as compare to F0 high sugar treatment in each sex?
    - what fraction of genes are expressed in opposite direction as compare to F0 high sugar treatment in each sex?
3. F2 generation:
    - What is the effect of F0 male parent diet(high sugar vs control) on the F2 progeny's gene expression pattern.
    - Is it F2 progeny sex-specific? 
    - what fraction of genes are expressed in same direction as compare to F0 high sugar treatment in each sex?
    - what fraction of genes are expressed in opposite direction as compare to F0 high sugar treatment in each sex?
4. If tryglcerides levels in all generation is in concordance with the gene expression profile? 

## For data visualization refer to this paper
* [limma powers differential expression analyses for RNA-sequencing and microarray studies](https://academic.oup.com/nar/article/43/7/e47/2414268)

## Key findings in EDA
### Number of genes which are differentially expressed (abs(FC) >= 1.3 and p-value < 0.05)
|Group       |  Total| up   |down  |
|----------  | -----: |-----:|-----:|
|f0_female   |  1845  | 856  |  989 |
|f0_male     |   728  | 249  |  479 |
|f1_femaleHS |   127  |  72  |   55 |
|f1_femaleLS |   206  | 181  |   25 |
|f1_maleHS   |    36  |   9  |   27 |
|f1_maleLS   |   181  | 128  |   53 |
|f2_femaleHS |   160  |  69  |   91 |
|f2_femaleLS |   242  |  84  |  158 |
|f2_maleHS   |   101  |  88  |   13 |
|f2_maleLS   |   484  | 407  |   77 |

## Note on GO enrichment analysis
1. Hypergeometric tests will be used for all enrichment analysis
2. The Gene ontology database or Affymetrix drosophila2.db database will be used as GO universe/background
3. Don't use DAVID since it's outdated
4. Decide which package should be used for analysis
    - ReactomePA/ Clusterprofiler
    - topGO
    - limma

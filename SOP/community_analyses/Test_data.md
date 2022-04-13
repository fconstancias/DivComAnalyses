---
title: ' Test dataset preparation 2022-04-12'
author: "Florentin CONSTANCIAS"
date: "2022-04-12"
output:
  html_document: 
    toc: yes
    toc_depth: 4
    keep_md: yes
---
<style type="text/css">
div.main-container {
max-width: 1800px;
margin-left: auto;
margin-right: auto;
}
.Table-Normal {
position: relative;
//display: block;
margin: 10px auto;
padding: 0;
width: 100%;
height: auto;
border-collapse: collapse;
text-align: center;
}
</style>




```r
require(tidyverse); packageVersion("tidyverse")
```

```
## Loading required package: tidyverse
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
## ✓ tibble  3.1.6     ✓ dplyr   1.0.7
## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
## ✓ readr   2.0.1     ✓ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
require(phyloseq); packageVersion("phyloseq")
```

```
## Loading required package: phyloseq
```

```r
options(getClass.msg=FALSE) # https://github.com/epurdom/clusterExperiment/issues/66
```

# PolyFermS:


```r
"~/Documents/GitHub/NRP72-FBT/data/processed/16S/ps_silva_dada2_human_chicken_meta.RDS" %>% 
  readRDS() %>% 
  subset_samples(Model == "Chicken" & Experiment == "Continuous" & Enrichment == "NotEnriched" & is.na(Paul) & Treatment %in% c("UNTREATED", "CTX", "VAN") & !Reactor %in% c("CR", "IR2")) %>%
  # subset_samples(tabled < 20000 & tabled > 1000) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_poly

ps_poly %>% 
  sample_data() %>% 
  data.frame() %>% 
  select(-I7_Index_ID, -index, -I5_Index_ID, -index2, -metagenomic_sample_name,-raw_metagenomic_pairs,-Reactor_Treatment_Dose ,-Treatment_Dose,-Model, -Date, -HV292.1_Copy_Number_permL, -Enrichment, -CCUG59168_Copy_Number_permL,-CTX_Copy_Number_permL, -VAN_Copy_Number_permL, -Reactor_Treatment, -Sample_description) %>% 
  mutate(Treatment = recode(Treatment, 
                            `UNTREATED` = "Control",
                            `CTX` = "Antibiotic_A",
                            `VAN` = "Antibiotic_B"),
         Antibiotic = recode(Antibiotic, 
                             `UNTREATED` = "Control",
                             `CTX` = "Antibiotic_A",
                             `VAN` = "Antibiotic_B")) -> sample_data(ps_poly)

ps_poly %>% 
  sample_data() %>% 
  data.frame()
```

```
##                  sample input filtered denoisedF denoisedR merged tabled
## IR1-16-S307 IR1-16-S307 16273    16228     16215     16126  15817  15817
## IR1-17-S279 IR1-17-S279 30887    30807     30773     30662  29870  29870
## IR1-18-S287 IR1-18-S287 35051    34926     34793     34693  33078  33078
## IR1-2-S96     IR1-2-S96 30159    30092     29944     29869  28243  28243
## IR1-24-S295 IR1-24-S295 14265    14243     14204     14153  13732  13732
## IR1-26-S119 IR1-26-S119 23432    23326     23307     23183  22620  22620
## IR1-3-S320   IR1-3-S320  8883     8865      8849      8809   8524   8524
## IR1-36-S159 IR1-36-S159 13340    13281     13264     13212  12937  12937
## IR1-37-S343 IR1-37-S343 21794    21690     21647     21601  21252  21252
## IR1-38-S337 IR1-38-S337 23007    22866     22847     22795  22423  22423
## IR1-39-S351 IR1-39-S351  8706     8655      8649      8623   8534   8534
## IR1-4-S223   IR1-4-S223 15037    15007     14984     14952  14505  14505
## IR1-43-S256 IR1-43-S256 17519    17412     17400     17329  16917  16917
## IR1-44-S264 IR1-44-S264 28866    28705     28687     28583  28096  28096
## IR1-45-S128 IR1-45-S128 11987    11921     11899     11862  11554  11554
## IR1-48-S339 IR1-48-S339 16625    16570     16550     16481  16111  16111
## IR1-54-S130 IR1-54-S130 26102    25956     25919     25869  25383  25383
## IR1-57-S356 IR1-57-S356 15767    15716     15688     15602  15202  15202
## IR1-6-S250   IR1-6-S250 29718    29662     29625     29553  28757  28757
## IR1-63-S155 IR1-63-S155 20842    20719     20690     20617  20157  20157
## IR1-69-S198 IR1-69-S198   132      131       128        91     72     72
## IR1-77-S202 IR1-77-S202 10448    10407     10391     10370  10153  10153
## IR1-8-S260   IR1-8-S260 25849    25766     25742     25634  25036  25036
## TR2-1-S273   TR2-1-S273 32834    32778     32729     32626  31606  31606
## TR2-13-S231 TR2-13-S231 22002    21915     21880     21840  21513  21513
## TR2-14-S232 TR2-14-S232 21133    21084     21028     20983  20461  20461
## TR2-15-S288 TR2-15-S288 31566    31486     31385     31307  30381  30381
## TR2-16-S184 TR2-16-S184 36760    36629     36533     36439  35449  35449
## TR2-17-S248 TR2-17-S248 31490    31404     31369     31205  30516  30516
## TR2-18-S326 TR2-18-S326 16339    16307     16275     16201  15845  15845
## TR2-19-S318 TR2-19-S318 10065    10054     10040     10015   9814   9814
## TR2-20-S271 TR2-20-S271 31072    31018     30982     30901  30167  30167
## TR2-21-S228 TR2-21-S228 15593    15567     15543     15456  14959  14959
## TR2-22-S364 TR2-22-S364 18760    18715     18688     18581  18015  18015
## TR2-25-S133 TR2-25-S133 10044    10028     10010      9980   9694   9694
## TR2-28-S286 TR2-28-S286 33858    33802     33626     33571  31718  31718
## TR2-3-S297   TR2-3-S297 16116    16088     16059     16017  15428  15428
## TR2-31-S182 TR2-31-S182 26488    26432     26391     26250  25319  25319
## TR2-34-S302 TR2-34-S302   114      113       108       100     84     84
## TR2-37-S330 TR2-37-S330 16482    16456     16434     16356  15811  15811
## TR2-40-S259 TR2-40-S259 21185    21122     21105     21025  20547  20547
## TR2-46-S239 TR2-46-S239 25929    25888     25862     25757  25098  25098
## TR2-52-S230 TR2-52-S230 14560    14526     14512     14459  14058  14058
## TR2-58-S172 TR2-58-S172 24670    24606     24565     24436  23412  23412
## TR2-64-S216 TR2-64-S216 21707    21663     21627     21526  20850  20850
## TR4-1-S193   TR4-1-S193 55187    54986     54963     54881  54034  54034
## TR4-13-S142 TR4-13-S142 19207    19074     19048     18973  18531  18531
## TR4-14-S226 TR4-14-S226 20906    20783     20752     20695  20340  20340
## TR4-15-S218 TR4-15-S218 21839    21709     21687     21638  21317  21317
## TR4-16-S276 TR4-16-S276 17118    17032     17006     16951  16595  16595
## TR4-17-S114 TR4-17-S114 28993    28752     28709     28624  28097  28097
## TR4-18-S200 TR4-18-S200   104      102        99        81     53     53
## TR4-19-S270 TR4-19-S270 28080    27827     27811     27758  27315  27315
## TR4-20-S269 TR4-20-S269 20518    20360     20353     20315  20002  20002
## TR4-21-S254 TR4-21-S254 21868    21735     21727     21684  21405  21405
## TR4-22-S344 TR4-22-S344 16291    16251     16227     16158  15906  15906
## TR4-25-S333 TR4-25-S333 14939    14883     14871     14817  14608  14608
## TR4-28-S183 TR4-28-S183 13668    13587     13562     13510  13345  13345
## TR4-3-S309   TR4-3-S309 22955    22915     22888     22810  22206  22206
## TR4-31-S278 TR4-31-S278 30109    29927     29902     29815  29549  29549
## TR4-34-S342 TR4-34-S342 13943    13880     13863     13793  13647  13647
## TR4-37-S188 TR4-37-S188 24208    24010     23977     23844  23524  23524
## TR4-40-S242 TR4-40-S242 33403    33159     33099     33044  32725  32725
## TR4-46-S215 TR4-46-S215 25294    25090     25053     24965  24729  24729
## TR4-52-S95   TR4-52-S95 35250    35030     34978     34883  34505  34505
## TR4-58-S112 TR4-58-S112 12291    12236     12219     12176  12075  12075
## TR4-64-S220 TR4-64-S220 15761    15686     15672     15607  15483  15483
##             filtered_pc denoisedF_pc denoisedR_pc merged_pc filtered_merged_pc
## IR1-16-S307       0.997        0.999        0.994     0.975              0.975
## IR1-17-S279       0.997        0.999        0.995     0.971              0.970
## IR1-18-S287       0.996        0.996        0.993     0.951              0.947
## IR1-2-S96         0.998        0.995        0.993     0.943              0.939
## IR1-24-S295       0.998        0.997        0.994     0.967              0.964
## IR1-26-S119       0.995        0.999        0.994     0.971              0.970
## IR1-3-S320        0.998        0.998        0.994     0.963              0.962
## IR1-36-S159       0.996        0.999        0.995     0.975              0.974
## IR1-37-S343       0.995        0.998        0.996     0.982              0.980
## IR1-38-S337       0.994        0.999        0.997     0.981              0.981
## IR1-39-S351       0.994        0.999        0.996     0.987              0.986
## IR1-4-S223        0.998        0.998        0.996     0.968              0.967
## IR1-43-S256       0.994        0.999        0.995     0.972              0.972
## IR1-44-S264       0.994        0.999        0.996     0.979              0.979
## IR1-45-S128       0.994        0.998        0.995     0.971              0.969
## IR1-48-S339       0.997        0.999        0.995     0.973              0.972
## IR1-54-S130       0.994        0.999        0.997     0.979              0.978
## IR1-57-S356       0.997        0.998        0.993     0.969              0.967
## IR1-6-S250        0.998        0.999        0.996     0.971              0.969
## IR1-63-S155       0.994        0.999        0.995     0.974              0.973
## IR1-69-S198       0.992        0.977        0.695     0.562              0.550
## IR1-77-S202       0.996        0.998        0.996     0.977              0.976
## IR1-8-S260        0.997        0.999        0.995     0.973              0.972
## TR2-1-S273        0.998        0.999        0.995     0.966              0.964
## TR2-13-S231       0.996        0.998        0.997     0.983              0.982
## TR2-14-S232       0.998        0.997        0.995     0.973              0.970
## TR2-15-S288       0.997        0.997        0.994     0.968              0.965
## TR2-16-S184       0.996        0.997        0.995     0.970              0.968
## TR2-17-S248       0.997        0.999        0.994     0.973              0.972
## TR2-18-S326       0.998        0.998        0.993     0.974              0.972
## TR2-19-S318       0.999        0.999        0.996     0.977              0.976
## TR2-20-S271       0.998        0.999        0.996     0.974              0.973
## TR2-21-S228       0.998        0.998        0.993     0.962              0.961
## TR2-22-S364       0.998        0.999        0.993     0.964              0.963
## TR2-25-S133       0.998        0.998        0.995     0.968              0.967
## TR2-28-S286       0.998        0.995        0.993     0.943              0.938
## TR2-3-S297        0.998        0.998        0.996     0.961              0.959
## TR2-31-S182       0.998        0.998        0.993     0.959              0.958
## TR2-34-S302       0.991        0.956        0.885     0.778              0.743
## TR2-37-S330       0.998        0.999        0.994     0.962              0.961
## TR2-40-S259       0.997        0.999        0.995     0.974              0.973
## TR2-46-S239       0.998        0.999        0.995     0.970              0.969
## TR2-52-S230       0.998        0.999        0.995     0.969              0.968
## TR2-58-S172       0.997        0.998        0.993     0.953              0.951
## TR2-64-S216       0.998        0.998        0.994     0.964              0.962
## TR4-1-S193        0.996        1.000        0.998     0.983              0.983
## TR4-13-S142       0.993        0.999        0.995     0.973              0.972
## TR4-14-S226       0.994        0.999        0.996     0.980              0.979
## TR4-15-S218       0.994        0.999        0.997     0.983              0.982
## TR4-16-S276       0.995        0.998        0.995     0.976              0.974
## TR4-17-S114       0.992        0.999        0.996     0.979              0.977
## TR4-18-S200       0.981        0.971        0.794     0.535              0.520
## TR4-19-S270       0.991        0.999        0.998     0.982              0.982
## TR4-20-S269       0.992        1.000        0.998     0.983              0.982
## TR4-21-S254       0.994        1.000        0.998     0.985              0.985
## TR4-22-S344       0.998        0.999        0.994     0.980              0.979
## TR4-25-S333       0.996        0.999        0.996     0.982              0.982
## TR4-28-S183       0.994        0.998        0.994     0.984              0.982
## TR4-3-S309        0.998        0.999        0.995     0.970              0.969
## TR4-31-S278       0.994        0.999        0.996     0.988              0.987
## TR4-34-S342       0.995        0.999        0.994     0.984              0.983
## TR4-37-S188       0.992        0.999        0.993     0.981              0.980
## TR4-40-S242       0.993        0.998        0.997     0.989              0.987
## TR4-46-S215       0.992        0.999        0.995     0.987              0.986
## TR4-52-S95        0.994        0.999        0.996     0.986              0.985
## TR4-58-S112       0.996        0.999        0.995     0.988              0.987
## TR4-64-S220       0.995        0.999        0.995     0.988              0.987
##             input_merged_pc tabled_joined chimera_out length_filtered tabled_pc
## IR1-16-S307           0.972         15817       15593           15593         1
## IR1-17-S279           0.967         29870       29056           29056         1
## IR1-18-S287           0.944         33078       29765           29765         1
## IR1-2-S96             0.936         28243       25344           25344         1
## IR1-24-S295           0.963         13732       13307           13307         1
## IR1-26-S119           0.965         22620       22248           22246         1
## IR1-3-S320            0.960          8524        8417            8417         1
## IR1-36-S159           0.970         12937       12724           12724         1
## IR1-37-S343           0.975         21252       20777           20772         1
## IR1-38-S337           0.975         22423       21618           21614         1
## IR1-39-S351           0.980          8534        8340            8339         1
## IR1-4-S223            0.965         14505       14293           14293         1
## IR1-43-S256           0.966         16917       16642           16642         1
## IR1-44-S264           0.973         28096       27507           27507         1
## IR1-45-S128           0.964         11554       11363           11363         1
## IR1-48-S339           0.969         16111       15811           15811         1
## IR1-54-S130           0.972         25383       24758           24758         1
## IR1-57-S356           0.964         15202       14958           14958         1
## IR1-6-S250            0.968         28757       28202           28202         1
## IR1-63-S155           0.967         20157       19681           19679         1
## IR1-69-S198           0.545            72          72              72         1
## IR1-77-S202           0.972         10153       10002           10000         1
## IR1-8-S260            0.969         25036       24358           24358         1
## TR2-1-S273            0.963         31606       30573           30573         1
## TR2-13-S231           0.978         21513       21064           21063         1
## TR2-14-S232           0.968         20461       19700           19700         1
## TR2-15-S288           0.962         30381       28832           28832         1
## TR2-16-S184           0.964         35449       34188           34185         1
## TR2-17-S248           0.969         30516       30071           30071         1
## TR2-18-S326           0.970         15845       15515           15515         1
## TR2-19-S318           0.975          9814        9662            9662         1
## TR2-20-S271           0.971         30167       29595           29595         1
## TR2-21-S228           0.959         14959       14630           14630         1
## TR2-22-S364           0.960         18015       17620           17620         1
## TR2-25-S133           0.965          9694        9527            9527         1
## TR2-28-S286           0.937         31718       28751           28751         1
## TR2-3-S297            0.957         15428       15094           15094         1
## TR2-31-S182           0.956         25319       24528           24528         1
## TR2-34-S302           0.737            84          84              84         1
## TR2-37-S330           0.959         15811       15473           15473         1
## TR2-40-S259           0.970         20547       20312           20312         1
## TR2-46-S239           0.968         25098       24475           24475         1
## TR2-52-S230           0.966         14058       13833           13833         1
## TR2-58-S172           0.949         23412       22828           22828         1
## TR2-64-S216           0.961         20850       20344           20344         1
## TR4-1-S193            0.979         54034       51873           51857         1
## TR4-13-S142           0.965         18531       18061           18061         1
## TR4-14-S226           0.973         20340       20007           20007         1
## TR4-15-S218           0.976         21317       20936           20936         1
## TR4-16-S276           0.969         16595       16308           16308         1
## TR4-17-S114           0.969         28097       27386           27381         1
## TR4-18-S200           0.510            53          52              52         1
## TR4-19-S270           0.973         27315       26523           26523         1
## TR4-20-S269           0.975         20002       19498           19494         1
## TR4-21-S254           0.979         21405       20811           20805         1
## TR4-22-S344           0.976         15906       15392           15392         1
## TR4-25-S333           0.978         14608       14213           14213         1
## TR4-28-S183           0.976         13345       12958           12958         1
## TR4-3-S309            0.967         22206       21693           21693         1
## TR4-31-S278           0.981         29549       28468           28465         1
## TR4-34-S342           0.979         13647       13296           13296         1
## TR4-37-S188           0.972         23524       22202           22197         1
## TR4-40-S242           0.980         32725       31714           31714         1
## TR4-46-S215           0.978         24729       23673           23673         1
## TR4-52-S95            0.979         34505       32412           32400         1
## TR4-58-S112           0.982         12075       11655           11655         1
## TR4-64-S220           0.982         15483       15160           15160         1
##             chimera_out_pc length_filtered_pc Description2 Experiment Reactor
## IR1-16-S307           0.99                  1         <NA> Continuous      IR
## IR1-17-S279           0.97                  1         <NA> Continuous      IR
## IR1-18-S287           0.90                  1         <NA> Continuous      IR
## IR1-2-S96             0.90                  1         <NA> Continuous      IR
## IR1-24-S295           0.97                  1         <NA> Continuous      IR
## IR1-26-S119           0.98                  1         <NA> Continuous      IR
## IR1-3-S320            0.99                  1         <NA> Continuous      IR
## IR1-36-S159           0.98                  1         <NA> Continuous      IR
## IR1-37-S343           0.98                  1         <NA> Continuous      IR
## IR1-38-S337           0.96                  1         <NA> Continuous      IR
## IR1-39-S351           0.98                  1         <NA> Continuous      IR
## IR1-4-S223            0.99                  1         <NA> Continuous      IR
## IR1-43-S256           0.98                  1         <NA> Continuous      IR
## IR1-44-S264           0.98                  1         <NA> Continuous      IR
## IR1-45-S128           0.98                  1         <NA> Continuous      IR
## IR1-48-S339           0.98                  1         <NA> Continuous      IR
## IR1-54-S130           0.98                  1         <NA> Continuous      IR
## IR1-57-S356           0.98                  1         <NA> Continuous      IR
## IR1-6-S250            0.98                  1         <NA> Continuous      IR
## IR1-63-S155           0.98                  1         <NA> Continuous      IR
## IR1-69-S198           1.00                  1         <NA> Continuous      IR
## IR1-77-S202           0.99                  1         <NA> Continuous      IR
## IR1-8-S260            0.97                  1         <NA> Continuous      IR
## TR2-1-S273            0.97                  1         <NA> Continuous     TR2
## TR2-13-S231           0.98                  1         <NA> Continuous     TR2
## TR2-14-S232           0.96                  1         <NA> Continuous     TR2
## TR2-15-S288           0.95                  1         <NA> Continuous     TR2
## TR2-16-S184           0.96                  1         <NA> Continuous     TR2
## TR2-17-S248           0.99                  1         <NA> Continuous     TR2
## TR2-18-S326           0.98                  1         <NA> Continuous     TR2
## TR2-19-S318           0.98                  1         <NA> Continuous     TR2
## TR2-20-S271           0.98                  1         <NA> Continuous     TR2
## TR2-21-S228           0.98                  1         <NA> Continuous     TR2
## TR2-22-S364           0.98                  1         <NA> Continuous     TR2
## TR2-25-S133           0.98                  1         <NA> Continuous     TR2
## TR2-28-S286           0.91                  1         <NA> Continuous     TR2
## TR2-3-S297            0.98                  1         <NA> Continuous     TR2
## TR2-31-S182           0.97                  1         <NA> Continuous     TR2
## TR2-34-S302           1.00                  1         <NA> Continuous     TR2
## TR2-37-S330           0.98                  1         <NA> Continuous     TR2
## TR2-40-S259           0.99                  1         <NA> Continuous     TR2
## TR2-46-S239           0.98                  1         <NA> Continuous     TR2
## TR2-52-S230           0.98                  1         <NA> Continuous     TR2
## TR2-58-S172           0.98                  1         <NA> Continuous     TR2
## TR2-64-S216           0.98                  1         <NA> Continuous     TR2
## TR4-1-S193            0.96                  1         <NA> Continuous     TR4
## TR4-13-S142           0.97                  1         <NA> Continuous     TR4
## TR4-14-S226           0.98                  1         <NA> Continuous     TR4
## TR4-15-S218           0.98                  1         <NA> Continuous     TR4
## TR4-16-S276           0.98                  1         <NA> Continuous     TR4
## TR4-17-S114           0.97                  1         <NA> Continuous     TR4
## TR4-18-S200           0.98                  1         <NA> Continuous     TR4
## TR4-19-S270           0.97                  1         <NA> Continuous     TR4
## TR4-20-S269           0.97                  1         <NA> Continuous     TR4
## TR4-21-S254           0.97                  1         <NA> Continuous     TR4
## TR4-22-S344           0.97                  1         <NA> Continuous     TR4
## TR4-25-S333           0.97                  1         <NA> Continuous     TR4
## TR4-28-S183           0.97                  1         <NA> Continuous     TR4
## TR4-3-S309            0.98                  1         <NA> Continuous     TR4
## TR4-31-S278           0.96                  1         <NA> Continuous     TR4
## TR4-34-S342           0.97                  1         <NA> Continuous     TR4
## TR4-37-S188           0.94                  1         <NA> Continuous     TR4
## TR4-40-S242           0.97                  1         <NA> Continuous     TR4
## TR4-46-S215           0.96                  1         <NA> Continuous     TR4
## TR4-52-S95            0.94                  1         <NA> Continuous     TR4
## TR4-58-S112           0.97                  1         <NA> Continuous     TR4
## TR4-64-S220           0.98                  1         <NA> Continuous     TR4
##                Treatment Day_of_Connection Day_of_Treatment Day_from_Inoculum
## IR1-16-S307      Control                -7              -23                16
## IR1-17-S279      Control                -6              -22                17
## IR1-18-S287      Control                -4              -20                19
## IR1-2-S96        Control               -21              -37                 2
## IR1-24-S295      Control                 1              -15                24
## IR1-26-S119      Control                 3              -13                26
## IR1-3-S320       Control               -20              -36                 3
## IR1-36-S159      Control                13               -3                36
## IR1-37-S343      Control                14               -2                37
## IR1-38-S337      Control                15               -1                38
## IR1-39-S351      Control                16                0                39
## IR1-4-S223       Control               -19              -35                 4
## IR1-43-S256      Control                20                4                43
## IR1-44-S264      Control                21                5                44
## IR1-45-S128      Control                22                6                45
## IR1-48-S339      Control                25                9                48
## IR1-54-S130      Control                31               19                54
## IR1-57-S356      Control                34               18                57
## IR1-6-S250       Control               -17              -33                 6
## IR1-63-S155      Control                40               24                63
## IR1-69-S198      Control                46               30                69
## IR1-77-S202      Control                54               38                77
## IR1-8-S260       Control               -15              -31                 8
## TR2-1-S273  Antibiotic_A                 1              -15                24
## TR2-13-S231 Antibiotic_A                13               -3                36
## TR2-14-S232 Antibiotic_A                14               -2                37
## TR2-15-S288 Antibiotic_A                15               -1                38
## TR2-16-S184 Antibiotic_A                16                0                39
## TR2-17-S248 Antibiotic_A                17                1                40
## TR2-18-S326 Antibiotic_A                18                2                41
## TR2-19-S318 Antibiotic_A                19                3                42
## TR2-20-S271 Antibiotic_A                20                4                43
## TR2-21-S228 Antibiotic_A                21                5                44
## TR2-22-S364 Antibiotic_A                22                6                45
## TR2-25-S133 Antibiotic_A                25                9                48
## TR2-28-S286 Antibiotic_A                28               12                51
## TR2-3-S297  Antibiotic_A                 3              -13                26
## TR2-31-S182 Antibiotic_A                31               15                54
## TR2-34-S302 Antibiotic_A                34               18                57
## TR2-37-S330 Antibiotic_A                37               21                60
## TR2-40-S259 Antibiotic_A                40               24                63
## TR2-46-S239 Antibiotic_A                46               30                69
## TR2-52-S230 Antibiotic_A                52               36                75
## TR2-58-S172 Antibiotic_A                58               42                81
## TR2-64-S216 Antibiotic_A                64               48                87
## TR4-1-S193  Antibiotic_B                 1              -15                24
## TR4-13-S142 Antibiotic_B                13               -3                36
## TR4-14-S226 Antibiotic_B                14               -2                37
## TR4-15-S218 Antibiotic_B                15               -1                38
## TR4-16-S276 Antibiotic_B                16                0                39
## TR4-17-S114 Antibiotic_B                17                1                40
## TR4-18-S200 Antibiotic_B                18                2                41
## TR4-19-S270 Antibiotic_B                19                3                42
## TR4-20-S269 Antibiotic_B                20                4                43
## TR4-21-S254 Antibiotic_B                21                5                44
## TR4-22-S344 Antibiotic_B                22                6                45
## TR4-25-S333 Antibiotic_B                25                9                48
## TR4-28-S183 Antibiotic_B                28               12                51
## TR4-3-S309  Antibiotic_B                 3              -13                26
## TR4-31-S278 Antibiotic_B                31               15                54
## TR4-34-S342 Antibiotic_B                34               18                57
## TR4-37-S188 Antibiotic_B                37               21                60
## TR4-40-S242 Antibiotic_B                40               24                63
## TR4-46-S215 Antibiotic_B                46               30                69
## TR4-52-S95  Antibiotic_B                52               36                75
## TR4-58-S112 Antibiotic_B                58               42                81
## TR4-64-S220 Antibiotic_B                64               48                87
##             Phase Treatment2 Paul GeneCopyNumberperML Antibiotic_mg.mL
## IR1-16-S307  Stab  UNTREATED <NA>            3.26e+10               NA
## IR1-17-S279  Stab  UNTREATED <NA>            3.20e+10               NA
## IR1-18-S287  Stab  UNTREATED <NA>            2.79e+10               NA
## IR1-2-S96    Stab  UNTREATED <NA>            2.60e+10               NA
## IR1-24-S295  Stab  UNTREATED <NA>            2.49e+10               NA
## IR1-26-S119  Stab  UNTREATED <NA>            2.82e+10               NA
## IR1-3-S320   Stab  UNTREATED <NA>            4.42e+10               NA
## IR1-36-S159  Stab  UNTREATED <NA>            3.15e+10               NA
## IR1-37-S343  Stab  UNTREATED <NA>            5.33e+10               NA
## IR1-38-S337  Stab  UNTREATED <NA>            6.54e+10               NA
## IR1-39-S351 Treat  UNTREATED <NA>            6.42e+10               NA
## IR1-4-S223   Stab  UNTREATED <NA>            4.44e+10               NA
## IR1-43-S256 Treat  UNTREATED <NA>            5.42e+10               NA
## IR1-44-S264 Treat  UNTREATED <NA>            4.10e+10               NA
## IR1-45-S128 Treat  UNTREATED <NA>            3.32e+10               NA
## IR1-48-S339 Treat  UNTREATED <NA>            2.95e+10               NA
## IR1-54-S130 Treat  UNTREATED <NA>            2.27e+10               NA
## IR1-57-S356 Treat  UNTREATED <NA>            2.04e+10               NA
## IR1-6-S250   Stab  UNTREATED <NA>            4.10e+10               NA
## IR1-63-S155 Treat  UNTREATED <NA>            2.56e+10               NA
## IR1-69-S198 Treat  UNTREATED <NA>            3.99e+10               NA
## IR1-77-S202 Treat  UNTREATED <NA>            3.14e+10               NA
## IR1-8-S260   Stab  UNTREATED <NA>            5.58e+10               NA
## TR2-1-S273   Stab         AB <NA>            3.49e+10               20
## TR2-13-S231  Stab         AB <NA>            5.27e+10               20
## TR2-14-S232  Stab         AB <NA>            5.78e+10               20
## TR2-15-S288  Stab         AB <NA>            6.23e+10               20
## TR2-16-S184 Treat         AB <NA>            4.77e+10               20
## TR2-17-S248 Treat         AB <NA>            2.71e+10               20
## TR2-18-S326 Treat         AB <NA>            3.22e+10               20
## TR2-19-S318 Treat         AB <NA>            4.20e+10               20
## TR2-20-S271 Treat         AB <NA>            2.27e+10               20
## TR2-21-S228 Treat         AB <NA>            2.89e+10               20
## TR2-22-S364 Treat         AB <NA>            3.98e+10               20
## TR2-25-S133 Treat         AB <NA>            3.16e+10               20
## TR2-28-S286 Treat         AB <NA>            3.16e+10               20
## TR2-3-S297   Stab         AB <NA>            4.23e+10               20
## TR2-31-S182 Treat         AB <NA>            2.61e+10               20
## TR2-34-S302 Treat         AB <NA>            2.79e+10               20
## TR2-37-S330 Treat         AB <NA>            2.08e+10               20
## TR2-40-S259 Treat         AB <NA>            2.77e+10               20
## TR2-46-S239 Treat         AB <NA>            5.05e+10               20
## TR2-52-S230 Treat         AB <NA>            3.52e+10               20
## TR2-58-S172 Treat         AB <NA>            2.54e+10               20
## TR2-64-S216 Treat         AB <NA>            2.79e+10               20
## TR4-1-S193   Stab         AB <NA>            3.40e+10               90
## TR4-13-S142  Stab         AB <NA>            5.63e+10               90
## TR4-14-S226  Stab         AB <NA>            5.00e+10               90
## TR4-15-S218  Stab         AB <NA>            6.68e+10               90
## TR4-16-S276 Treat         AB <NA>            3.96e+10               90
## TR4-17-S114 Treat         AB <NA>            3.74e+10               90
## TR4-18-S200 Treat         AB <NA>            6.03e+10               90
## TR4-19-S270 Treat         AB <NA>            5.61e+10               90
## TR4-20-S269 Treat         AB <NA>            4.27e+10               90
## TR4-21-S254 Treat         AB <NA>            5.82e+10               90
## TR4-22-S344 Treat         AB <NA>            7.52e+10               90
## TR4-25-S333 Treat         AB <NA>            4.16e+10               90
## TR4-28-S183 Treat         AB <NA>            4.30e+10               90
## TR4-3-S309   Stab         AB <NA>            3.98e+10               90
## TR4-31-S278 Treat         AB <NA>            4.38e+10               90
## TR4-34-S342 Treat         AB <NA>            4.68e+10               90
## TR4-37-S188 Treat         AB <NA>            5.12e+10               90
## TR4-40-S242 Treat         AB <NA>            5.57e+10               90
## TR4-46-S215 Treat         AB <NA>            6.31e+10               90
## TR4-52-S95  Treat         AB <NA>            4.94e+10               90
## TR4-58-S112 Treat         AB <NA>            7.74e+10               90
## TR4-64-S220 Treat         AB <NA>            6.94e+10               90
##             Fermentation   Antibiotic Lactose_mM Glucose_mM Galactose_mM
## IR1-16-S307           NA         <NA>      0.111      0.000        0.000
## IR1-17-S279           NA         <NA>      0.108      0.000        0.000
## IR1-18-S287           NA         <NA>      0.112      0.000        0.000
## IR1-2-S96             NA         <NA>      0.000      0.000        0.000
## IR1-24-S295           NA         <NA>      0.000      0.000        0.000
## IR1-26-S119           NA         <NA>      0.000      0.000        0.000
## IR1-3-S320            NA         <NA>      0.000      0.000        0.000
## IR1-36-S159           NA         <NA>      0.000      0.000        0.000
## IR1-37-S343           NA         <NA>      0.000         NA        0.000
## IR1-38-S337           NA         <NA>      0.000      0.000        1.608
## IR1-39-S351           NA         <NA>      0.000      0.000        1.765
## IR1-4-S223            NA         <NA>      0.000      0.736        0.000
## IR1-43-S256           NA         <NA>      0.000      0.000        0.000
## IR1-44-S264           NA         <NA>      0.000      0.000        0.000
## IR1-45-S128           NA         <NA>         NA         NA           NA
## IR1-48-S339           NA         <NA>      0.000      0.000        0.000
## IR1-54-S130           NA         <NA>      0.000      0.000        0.000
## IR1-57-S356           NA         <NA>      0.000      0.000        0.000
## IR1-6-S250            NA         <NA>      0.000      0.000        0.000
## IR1-63-S155           NA         <NA>      0.000      0.000        0.000
## IR1-69-S198           NA         <NA>      0.000      0.000        0.000
## IR1-77-S202           NA         <NA>      0.000      0.000        0.000
## IR1-8-S260            NA         <NA>      0.000      0.000        0.000
## TR2-1-S273            NA Antibiotic_A         NA         NA           NA
## TR2-13-S231           NA Antibiotic_A      0.000      0.000        0.000
## TR2-14-S232           NA Antibiotic_A      0.000      0.000        0.000
## TR2-15-S288           NA Antibiotic_A      0.000      0.000        0.000
## TR2-16-S184           NA Antibiotic_A      0.000      0.000        0.000
## TR2-17-S248           NA Antibiotic_A      0.000      0.000        0.000
## TR2-18-S326           NA Antibiotic_A      0.000      0.000        0.000
## TR2-19-S318           NA Antibiotic_A      0.000      0.000        0.000
## TR2-20-S271           NA Antibiotic_A      0.000      0.000        0.000
## TR2-21-S228           NA Antibiotic_A      0.000      0.000        0.000
## TR2-22-S364           NA Antibiotic_A      0.000      0.000        0.000
## TR2-25-S133           NA Antibiotic_A      0.000      0.000        0.000
## TR2-28-S286           NA Antibiotic_A      0.000      0.000        0.000
## TR2-3-S297            NA Antibiotic_A      0.000      0.000        0.000
## TR2-31-S182           NA Antibiotic_A      0.000      0.000        1.008
## TR2-34-S302           NA Antibiotic_A      0.000      0.000        0.000
## TR2-37-S330           NA Antibiotic_A      0.000      0.000        0.000
## TR2-40-S259           NA Antibiotic_A      0.000      0.000        0.000
## TR2-46-S239           NA Antibiotic_A      0.000      0.000        0.000
## TR2-52-S230           NA Antibiotic_A      0.000      0.000        0.000
## TR2-58-S172           NA Antibiotic_A      0.000      0.000        0.000
## TR2-64-S216           NA Antibiotic_A      0.000      0.000        0.000
## TR4-1-S193            NA Antibiotic_B         NA         NA           NA
## TR4-13-S142           NA Antibiotic_B      0.000      0.000        0.000
## TR4-14-S226           NA Antibiotic_B      0.000      0.000        0.000
## TR4-15-S218           NA Antibiotic_B      0.000      0.000        0.000
## TR4-16-S276           NA Antibiotic_B      0.000      0.000        0.000
## TR4-17-S114           NA Antibiotic_B      0.000      1.401        1.045
## TR4-18-S200           NA Antibiotic_B      0.000      0.000        1.532
## TR4-19-S270           NA Antibiotic_B      0.000      0.000        1.380
## TR4-20-S269           NA Antibiotic_B      0.000      0.446        1.237
## TR4-21-S254           NA Antibiotic_B      0.598      0.000        1.033
## TR4-22-S344           NA Antibiotic_B      0.471      0.383        0.000
## TR4-25-S333           NA Antibiotic_B      0.318      0.000        0.000
## TR4-28-S183           NA Antibiotic_B      0.000      0.000        0.000
## TR4-3-S309            NA Antibiotic_B      0.000      0.000        0.000
## TR4-31-S278           NA Antibiotic_B         NA         NA           NA
## TR4-34-S342           NA Antibiotic_B      0.000      0.000        0.000
## TR4-37-S188           NA Antibiotic_B      0.000      0.000        0.000
## TR4-40-S242           NA Antibiotic_B      0.000      0.000        0.000
## TR4-46-S215           NA Antibiotic_B      0.000      0.000        0.000
## TR4-52-S95            NA Antibiotic_B      0.000      0.000        0.000
## TR4-58-S112           NA Antibiotic_B      0.000      0.000        0.000
## TR4-64-S220           NA Antibiotic_B      0.000      0.000        0.000
##             Succinat_mM Lactat_mM Formiat_mM Acetat_mM Propionat_mM
## IR1-16-S307       8.403     0.000      3.022    59.445       14.806
## IR1-17-S279       6.582     0.000      2.551    63.570       12.533
## IR1-18-S287       6.977     0.000      3.131    61.230       12.129
## IR1-2-S96        11.188     0.000      2.195    54.634        4.767
## IR1-24-S295       6.862     0.000      3.239    65.522       15.644
## IR1-26-S119       4.957     0.000      0.000    70.004       16.081
## IR1-3-S320       12.546     0.000      0.000    65.405       10.260
## IR1-36-S159       1.818     0.000      0.000    80.582       13.261
## IR1-37-S343       2.801     7.732      0.000    50.595       10.602
## IR1-38-S337       3.233    10.316      3.202    45.582        9.269
## IR1-39-S351       4.742    11.710      3.238    40.536        7.296
## IR1-4-S223       13.888     0.000      0.000    74.813       12.634
## IR1-43-S256       0.696     0.000      0.000    84.344       18.354
## IR1-44-S264       2.068     0.000      0.000    76.906       16.969
## IR1-45-S128          NA        NA         NA        NA           NA
## IR1-48-S339       0.674     0.000      0.000    69.611       10.977
## IR1-54-S130       0.646     0.000      0.000    65.814       11.962
## IR1-57-S356       0.330     0.000      0.000    67.957       10.292
## IR1-6-S250       14.292     0.000      0.000    83.221       12.127
## IR1-63-S155       0.356     0.000      0.000    71.316       12.482
## IR1-69-S198       0.000     0.000      0.000    78.578       12.904
## IR1-77-S202       0.000     4.346      0.000    74.352       13.977
## IR1-8-S260       11.229     0.000      0.000    77.668       12.081
## TR2-1-S273           NA        NA         NA        NA           NA
## TR2-13-S231       1.990     0.000      0.000    93.560       14.283
## TR2-14-S232       1.687     0.000      0.000    90.856       13.870
## TR2-15-S288       0.792     0.000      0.000    98.070       13.180
## TR2-16-S184       0.752     0.000      0.000    99.902       12.294
## TR2-17-S248       1.116     0.000      0.000    95.235       12.018
## TR2-18-S326       2.193     0.000      0.000    96.527        9.384
## TR2-19-S318       2.973     0.000      0.000    90.991        7.423
## TR2-20-S271       3.025     0.000      0.000    75.634        7.726
## TR2-21-S228       3.420     0.000      0.000    84.906        5.716
## TR2-22-S364       2.853     0.000      0.000    84.486        4.727
## TR2-25-S133       3.076     0.000      0.000    82.344        5.823
## TR2-28-S286       0.877     0.000      0.000    88.906        4.666
## TR2-3-S297        6.593     0.000      2.587    63.021       18.800
## TR2-31-S182       2.177     0.000      0.000    91.660        8.217
## TR2-34-S302       4.281     0.000      0.000    94.764        5.288
## TR2-37-S330       3.602     0.000      0.000   103.752        4.216
## TR2-40-S259       5.119     0.000      0.000   102.943        8.903
## TR2-46-S239       3.776     0.000      0.000   110.663       12.096
## TR2-52-S230       0.000     0.000      0.000    91.397       12.261
## TR2-58-S172       0.327     0.000      1.249    68.746       13.115
## TR2-64-S216       0.511     0.000      2.219    57.910       13.311
## TR4-1-S193           NA        NA         NA        NA           NA
## TR4-13-S142       0.623     0.000      0.000    89.036       13.196
## TR4-14-S226       0.748     0.000      0.000    94.132       13.831
## TR4-15-S218       0.621     0.000      0.000   101.441       14.850
## TR4-16-S276       0.742     0.000      0.000    95.815       13.114
## TR4-17-S114       4.792     0.000      0.000    79.243       17.005
## TR4-18-S200       8.744     0.000      0.000    57.244       21.131
## TR4-19-S270      11.605     0.000      0.000    48.991       21.902
## TR4-20-S269      10.699     0.000     10.349    46.502       23.092
## TR4-21-S254      10.796     0.000      9.618    47.637       24.239
## TR4-22-S344      10.667     0.000      0.000    44.945       22.242
## TR4-25-S333       8.354     0.000      0.000    39.160       17.485
## TR4-28-S183       7.903     0.000      0.000    37.823       14.292
## TR4-3-S309        4.884     0.000      0.000    74.419       16.121
## TR4-31-S278          NA        NA         NA        NA           NA
## TR4-34-S342       5.187     0.000      0.000    44.793       15.056
## TR4-37-S188       6.401     0.000      0.000    44.769       16.874
## TR4-40-S242       7.825     7.755      0.000    42.346       15.733
## TR4-46-S215       7.444     8.524      0.000    46.362       19.358
## TR4-52-S95        7.698     0.000      0.000    40.062       23.478
## TR4-58-S112       6.615     0.000      0.000    30.977       15.826
## TR4-64-S220       0.418     3.808      0.000    41.394       17.022
##             Isobutyrat_mM Butyrat_mM Isovalerat_mM Valerat_mM Total_SCFA_mM
## IR1-16-S307         6.315     35.967         7.097      6.126       141.181
## IR1-17-S279         6.502     42.810         6.828      5.438       146.814
## IR1-18-S287         6.518     41.377         7.063      5.312       143.737
## IR1-2-S96           1.011     32.783         0.986      5.139       112.703
## IR1-24-S295         6.740     43.375         7.831      5.631       154.844
## IR1-26-S119         7.425     42.255         7.986      7.521       156.229
## IR1-3-S320          5.111     39.518         7.536      6.546       146.922
## IR1-36-S159         6.892     41.525         6.992      8.255       159.325
## IR1-37-S343         4.273     20.076         3.437      2.288       101.804
## IR1-38-S337         4.210     14.028         2.306      1.122        94.876
## IR1-39-S351         0.000      9.521         1.415      0.477        80.700
## IR1-4-S223          6.540     37.590         8.385      7.366       161.216
## IR1-43-S256         7.196     27.685         4.619      0.000       142.894
## IR1-44-S264         6.771     27.800         4.159      1.947       136.620
## IR1-45-S128            NA         NA            NA         NA            NA
## IR1-48-S339         9.553     47.429         4.358      6.721       149.323
## IR1-54-S130         9.046     41.833         5.919      6.912       142.132
## IR1-57-S356         5.339     39.561         4.819      6.555       134.853
## IR1-6-S250          6.863     34.543         7.530      6.715       165.291
## IR1-63-S155         6.156     39.345         5.694      6.666       142.015
## IR1-69-S198         6.546     40.828         5.852      8.191       152.899
## IR1-77-S202        10.539     40.389         3.007      7.588       154.198
## IR1-8-S260          6.830     31.547         7.009      4.427       150.791
## TR2-1-S273             NA         NA            NA         NA            NA
## TR2-13-S231         7.367     41.719         6.875      6.376       172.170
## TR2-14-S232         7.782     43.157         6.980      7.847       172.179
## TR2-15-S288         8.776     44.333         6.422      7.313       178.886
## TR2-16-S184         7.596     44.129         5.802      6.606       177.081
## TR2-17-S248         7.794     47.890         6.038      3.410       173.501
## TR2-18-S326         6.547     47.253         3.761      1.469       167.134
## TR2-19-S318         0.000     44.651         2.036      0.036       148.110
## TR2-20-S271         0.000     48.862         3.874      0.570       139.691
## TR2-21-S228         0.000     41.240         1.044      0.000       136.326
## TR2-22-S364         0.000     43.389         1.149      0.000       136.604
## TR2-25-S133         0.000     40.410         1.342      0.000       132.995
## TR2-28-S286         0.000     44.515         1.109      0.000       140.073
## TR2-3-S297          6.150     45.191         7.968      3.384       153.694
## TR2-31-S182         5.115     45.407         4.247      0.000       156.823
## TR2-34-S302         0.000     39.231         1.511      0.000       145.075
## TR2-37-S330         1.411     39.604         1.355      0.000       153.940
## TR2-40-S259         0.000     35.693         4.456      0.000       157.114
## TR2-46-S239         5.428     44.934         6.609      0.000       183.506
## TR2-52-S230         8.766     43.145         5.232      0.000       160.801
## TR2-58-S172        10.398     52.683         7.211      0.338       154.067
## TR2-64-S216         7.440     35.756         6.925      0.233       124.305
## TR4-1-S193             NA         NA            NA         NA            NA
## TR4-13-S142         7.595     36.478         6.940      8.164       162.032
## TR4-14-S226         8.933     34.897         6.840      8.443       167.824
## TR4-15-S218         7.813     33.872         6.542      8.752       173.891
## TR4-16-S276         8.162     37.582         6.557      8.358       170.330
## TR4-17-S114         8.207     28.479         6.182      4.307       148.215
## TR4-18-S200         5.169     14.109         5.840      1.585       113.822
## TR4-19-S270         4.715      6.737         5.970      0.436       100.356
## TR4-20-S269         4.802      3.750         5.741      0.000       104.935
## TR4-21-S254         0.000      1.522         5.250      0.000        99.062
## TR4-22-S344         6.994      2.367         5.154      0.000        92.369
## TR4-25-S333         5.935      2.601         5.816      0.000        79.351
## TR4-28-S183         5.811      2.036         4.958      0.000        72.823
## TR4-3-S309          6.230     45.385         8.218      3.351       158.608
## TR4-31-S278            NA         NA            NA         NA            NA
## TR4-34-S342         4.904      0.634         4.836      0.000        75.410
## TR4-37-S188         4.021      0.450         6.832      0.000        79.347
## TR4-40-S242         4.230      0.000         6.867      0.000        77.001
## TR4-46-S215         4.941      0.287         7.161      0.000        85.553
## TR4-52-S95         12.070      0.000         3.823      0.000        87.131
## TR4-58-S112         0.000      0.994         5.581      0.000        59.993
## TR4-64-S220         5.820      2.230         4.591      0.000        71.475
```


```r
ps_poly %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  saveRDS("~/Documents/GitHub/DivComAnalyses/data-raw/ps_PolyFermS.RDS")
```


# Invivo:

```r
"~/Documents/GitHub/mIMT/data/physeq_silva132_clean_phangorn_tree_NAin.RDS" %>% 
  readRDS() %>% 
  subset_samples(sample_type %in% c("Feces","Cecum")) %>%
  subset_samples(!treatment %in% c("InoculumFermenter_slurry","InoculumFermenter", "Caecal_inoculum_F6")) %>%
  subset_samples(!is.na(treatment)) %>%
  subset_samples(!is.na(day)) %>%
  subset_samples(remove != "x") %>% #sample_data() %>% data.frame()
  subset_samples(tabled > 500) %>% 
  filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> physeq
```


```r
sample_data(physeq)$treatment = factor(sample_data(physeq)$treatment,
                                       levels = c("H2O", "H2O_co",
                                                  "Eubiotic",
                                                  "Normobiotic",
                                                  "Dysbiotic",
                                                  "DSS", "DSS_co",
                                                  "Predni"))

sample_data(physeq)$sample_type = factor(sample_data(physeq)$sample_type,
                                         levels = c("Feces",
                                                    "Cecum"))

sample_data(physeq)$day = factor(sample_data(physeq)$day,
                                 levels = c("Day_-1",
                                            "Day_6",
                                            "Day_12",
                                            "Day_14"))

sample_data(physeq)$treatment_grouped <- factor(ifelse(sample_data(physeq)$treatment == "DSS_co", "DSS", 
                                                       ifelse(sample_data(physeq)$treatment == "H2O_co", "H2O", 
                                                              as.vector(sample_data(physeq)$treatment))),
                                                levels = c("H2O",
                                                           "Eubiotic",
                                                           "Normobiotic",
                                                           "Dysbiotic",
                                                           "DSS",
                                                           "Predni"))


sample_data(physeq)$treatment_invivo <- factor(ifelse(sample_data(physeq)$treatment_grouped %in% c("H2O") | sample_data(physeq)$day %in% c("Day_-1"), 
                                                      "none", 
                                                      as.vector(sample_data(physeq)$treatment_grouped)),
                                               levels = c("none",
                                                          "Eubiotic",
                                                          "Normobiotic",
                                                          "Dysbiotic",
                                                          "DSS",
                                                          "Predni"))

sample_data(physeq)$day_num <- as.numeric(parse_number(as.character(sample_data(physeq)$day)))

tax_table(physeq) <- tax_table(physeq)[,-8]


physeq %>% 
  sample_data() %>% 
  data.frame()
```

```
##             sample_file  input filtered filtered_pc denoisedF denoisedR
## MSQ4_AG_012 MSQ4_AG_012  93567    93211        1.00     92469     92315
## MSQ4_AG_013 MSQ4_AG_013  40818    40666        1.00     40327     40227
## MSQ4_AG_014 MSQ4_AG_014  53457    53157        0.99     52395     52446
## MSQ4_AG_015 MSQ4_AG_015  55794    55469        0.99     54765     54579
## MSQ4_AG_016 MSQ4_AG_016  40971    40388        0.99     39465     39648
## MSQ4_AG_017 MSQ4_AG_017  52401    52106        0.99     51510     51336
## MSQ4_AG_018 MSQ4_AG_018  49906    49796        1.00     49291     49546
## MSQ4_AG_019 MSQ4_AG_019  50731    50522        1.00     50000     49963
## MSQ4_AG_020 MSQ4_AG_020  69677    69339        1.00     68796     68640
## MSQ4_AG_021 MSQ4_AG_021  34990    34902        1.00     34611     34636
## MSQ4_AG_022 MSQ4_AG_022  55591    55355        1.00     54835     54735
## MSQ4_AG_023 MSQ4_AG_023  53398    53114        0.99     52558     52451
## MSQ4_AG_024 MSQ4_AG_024  46684    46414        0.99     45811     45946
## MSQ4_AG_025 MSQ4_AG_025  58371    57920        0.99     57113     57230
## MSQ4_AG_027 MSQ4_AG_027  59798    59534        1.00     58847     58892
## MSQ4_AG_028 MSQ4_AG_028  44799    44424        0.99     43891     43482
## MSQ4_AG_029 MSQ4_AG_029  39558    39318        0.99     38996     38835
## MSQ4_AG_030 MSQ4_AG_030  75323    75034        1.00     73911     74263
## MSQ4_AG_031 MSQ4_AG_031  57329    56983        0.99     56039     55936
## MSQ4_AG_033 MSQ4_AG_033  72638    71967        0.99     71090     71253
## MSQ4_AG_034 MSQ4_AG_034  53259    53048        1.00     52446     52420
## MSQ4_AG_035 MSQ4_AG_035  65565    65270        1.00     64634     64650
## MSQ4_AG_036 MSQ4_AG_036  77111    76640        0.99     75958     75085
## MSQ4_AG_037 MSQ4_AG_037  36922    36800        1.00     36357     36452
## MSQ4_AG_038 MSQ4_AG_038  53327    53112        1.00     52379     52522
## MSQ4_AG_039 MSQ4_AG_039  48920    48737        1.00     47846     48215
## MSQ4_AG_040 MSQ4_AG_040  38760    38291        0.99     37527     37748
## MSQ4_AG_041 MSQ4_AG_041  45415    44890        0.99     44161     44437
## MSQ4_AG_042 MSQ4_AG_042  52962    52710        1.00     52121     52200
## MSQ4_AG_043 MSQ4_AG_043  52039    51598        0.99     50997     51089
## MSQ4_AG_044 MSQ4_AG_044  52547    52402        1.00     51864     52036
## MSQ4_AG_045 MSQ4_AG_045  38678    38540        1.00     38170     38178
## MSQ4_AG_046 MSQ4_AG_046  50942    50761        1.00     49987     50055
## MSQ4_AG_047 MSQ4_AG_047  51012    50547        0.99     49780     49609
## MSQ4_AG_049 MSQ4_AG_049  60539    59726        0.99     58269     58706
## MSQ4_AG_050 MSQ4_AG_050  52454    52136        0.99     51273     50886
## MSQ4_AG_051 MSQ4_AG_051  58408    58062        0.99     56836     57236
## MSQ4_AG_052 MSQ4_AG_052  63806    63282        0.99     62447     62392
## MSQ4_AG_053 MSQ4_AG_053  45783    45404        0.99     44683     44819
## MSQ4_AG_054 MSQ4_AG_054  62902    62235        0.99     60999     61371
## MSQ4_AG_055 MSQ4_AG_055  62301    61916        0.99     60471     61055
## MSQ4_AG_056 MSQ4_AG_056  35572    34288        0.96     32717     33719
## MSQ4_AG_057 MSQ4_AG_057  63920    63563        0.99     62662     62913
## MSQ4_AG_058 MSQ4_AG_058  60572    60229        0.99     59580     59050
## MSQ4_AG_059 MSQ4_AG_059  61027    60756        1.00     60154     60162
## MSQ4_AG_062 MSQ4_AG_062  91051    90544        0.99     89611     88397
## MSQ4_AG_064 MSQ4_AG_064  49059    48696        0.99     47885     47755
## MSQ4_AG_065 MSQ4_AG_065  50683    50339        0.99     49570     49552
## MSQ4_AG_072 MSQ4_AG_072  27975    27845        1.00     27527     27312
## MSQ4_AG_076 MSQ4_AG_076  62101    61803        1.00     61320     60640
## MSQ4_AG_077 MSQ4_AG_077  48039    47769        0.99     47404     45592
## MSQ4_AG_080 MSQ4_AG_080  42443    41953        0.99     41150     41277
## MSQ4_AG_083 MSQ4_AG_083  55840    55409        0.99     55006     54267
## MSQ4_AG_084 MSQ4_AG_084  64132    63673        0.99     63328     62957
## MSQ4_AG_085 MSQ4_AG_085  48367    48121        0.99     47654     47338
## MSQ4_AG_086 MSQ4_AG_086  69480    69020        0.99     68163     67904
## MSQ4_AG_087 MSQ4_AG_087  50207    49895        0.99     49057     49414
## MSQ4_AG_088 MSQ4_AG_088  40287    39949        0.99     39410     39736
## MSQ4_AG_089 MSQ4_AG_089  59156    58962        1.00     58378     58546
## MSQ4_AG_091 MSQ4_AG_091  59808    59573        1.00     59038     59223
## MSQ4_AG_092 MSQ4_AG_092  66878    66665        1.00     66329     66184
## MSQ4_AG_093 MSQ4_AG_093  49707    49503        1.00     49236     49074
## MSQ4_AG_095 MSQ4_AG_095  53637    53374        1.00     52612     52473
## MSQ4_AG_096 MSQ4_AG_096  58366    58124        1.00     57558     57557
## MSQ4_AG_101 MSQ4_AG_101  44687    44476        1.00     43991     43442
## MSQ4_AG_103 MSQ4_AG_103  43368    43164        1.00     42566     42450
## MSQ4_AG_104 MSQ4_AG_104  42000    41770        0.99     41368     41071
## MSQ4_AG_105 MSQ4_AG_105  44824    44673        1.00     44195     43993
## MSQ4_AG_106 MSQ4_AG_106  27532    27427        1.00     27140     26935
## MSQ4_AG_107 MSQ4_AG_107  42230    42049        1.00     41694     41013
## MSQ4_AG_108 MSQ4_AG_108  39332    39034        0.99     38634     38446
## MSQ4_AG_109 MSQ4_AG_109  42091    41872        0.99     41455     41088
## MSQ4_AG_111 MSQ4_AG_111  40620    40433        1.00     40081     39347
## MSQ4_AG_112 MSQ4_AG_112  53189    53039        1.00     52770     52559
## MSQ4_AG_113 MSQ4_AG_113  47890    47651        1.00     47206     46771
## MSQ4_AG_117 MSQ4_AG_117  47191    46890        0.99     46426     45967
## MSQ4_AG_119 MSQ4_AG_119  40226    40027        1.00     39649     38938
## MSQ4_AG_121 MSQ4_AG_121  55552    55098        0.99     54488     54043
## MSQ4_AG_125 MSQ4_AG_125  39530    39256        0.99     38814     38782
## MSQ4_AG_127 MSQ4_AG_127  40641    40473        1.00     39890     39801
## MSQ4_AG_128 MSQ4_AG_128  62334    62031        1.00     61579     60901
## MSQ4_AG_129 MSQ4_AG_129  51425    51231        1.00     50801     50419
## MSQ4_AG_130 MSQ4_AG_130  44513    44300        1.00     43931     43765
## MSQ4_AG_131 MSQ4_AG_131  46907    46588        0.99     46174     45791
## MSQ4_AG_132 MSQ4_AG_132  62755    62478        1.00     61813     61613
## MSQ4_AG_133 MSQ4_AG_133  52811    52750        1.00     52608     52582
## MSQ4_AG_135 MSQ4_AG_135  48184    47986        1.00     47510     47650
## MSQ4_AG_136 MSQ4_AG_136  50255    50062        1.00     49589     49481
## MSQ4_AG_141 MSQ4_AG_141  42603    42338        0.99     41868     41462
## MSQ4_AG_143 MSQ4_AG_143  51996    51719        0.99     51228     51116
## MSQ4_AG_144 MSQ4_AG_144  56545    56149        0.99     55566     55094
## MSQ4_AG_145 MSQ4_AG_145  62003    61573        0.99     60862     59759
## MSQ4_AG_148 MSQ4_AG_148 106414   105300        0.99    103679    103442
## MSQ4_AG_149 MSQ4_AG_149  65720    65374        0.99     64492     64683
## MSQ4_AG_151 MSQ4_AG_151  51698    51461        1.00     50592     50846
## MSQ4_AG_152 MSQ4_AG_152  60137    59777        0.99     59149     58604
## MSQ4_AG_153 MSQ4_AG_153  74050    73761        1.00     73032     72948
## MSQ4_AG_154 MSQ4_AG_154  59604    59382        1.00     58712     58688
## MSQ4_AG_155 MSQ4_AG_155  51161    50935        1.00     50293     50069
## MSQ4_AG_156 MSQ4_AG_156  51215    50963        1.00     50328     50259
## MSQ4_AG_157 MSQ4_AG_157  52562    52278        0.99     51669     51628
## MSQ4_AG_159 MSQ4_AG_159  50577    50372        1.00     49820     49927
## MSQ4_AG_160 MSQ4_AG_160  52036    51881        1.00     51277     51128
## MSQ4_AG_161 MSQ4_AG_161  53262    53066        1.00     52382     52295
## MSQ4_AG_162 MSQ4_AG_162  26849    26741        1.00     26537     26404
## MSQ4_AG_163 MSQ4_AG_163  41089    40983        1.00     40544     40527
## MSQ4_AG_164 MSQ4_AG_164  58504    58260        1.00     57600     57629
## MSQ4_AG_165 MSQ4_AG_165  44066    43788        0.99     43433     43258
## MSQ4_AG_167 MSQ4_AG_167  40666    40526        1.00     40127     39915
## MSQ4_AG_168 MSQ4_AG_168  59813    59641        1.00     58927     59222
## MSQ4_AG_169 MSQ4_AG_169  52568    52350        1.00     51774     51678
## MSQ4_AG_170 MSQ4_AG_170  53081    52864        1.00     52362     52188
## MSQ4_AG_171 MSQ4_AG_171  54893    54774        1.00     54406     54406
## MSQ4_AG_172 MSQ4_AG_172  50625    50407        1.00     49826     49943
## MSQ4_AG_173 MSQ4_AG_173  53181    53027        1.00     52527     52632
## MSQ4_AG_175 MSQ4_AG_175  64177    63956        1.00     63311     63418
## MSQ4_AG_176 MSQ4_AG_176  52979    52633        0.99     52230     51930
## MSQ4_AG_177 MSQ4_AG_177  45011    44786        1.00     44214     44211
## MSQ4_AG_178 MSQ4_AG_178  39193    39021        1.00     38638     38559
## MSQ4_AG_179 MSQ4_AG_179  40775    40559        0.99     40207     39971
## MSQ4_AG_180 MSQ4_AG_180  42071    41887        1.00     41491     41500
## MSQ4_AG_181 MSQ4_AG_181  41096    40910        1.00     40546     40383
## MSQ4_AG_183 MSQ4_AG_183  40995    40747        0.99     40340     40184
## MSQ4_AG_184 MSQ4_AG_184  49569    49326        1.00     48760     48775
## MSQ4_AG_185 MSQ4_AG_185  48720    48497        1.00     47958     47940
## MSQ4_AG_186 MSQ4_AG_186  47408    47208        1.00     46771     46814
## MSQ4_AG_187 MSQ4_AG_187  37512    37373        1.00     36828     36760
## MSQ4_AG_188 MSQ4_AG_188  39871    39566        0.99     39104     39056
## MSQ4_AG_189 MSQ4_AG_189  56051    55723        0.99     55092     54971
## MSQ4_AG_191 MSQ4_AG_191  59032    58519        0.99     57455     57796
## MSQ4_AG_192 MSQ4_AG_192  46930    46753        1.00     46431     46398
## MSQ4_AG_193 MSQ4_AG_193  47507    47344        1.00     47037     46987
## MSQ4_AG_194 MSQ4_AG_194  44699    44583        1.00     44407     43977
## MSQ4_AG_195 MSQ4_AG_195  36586    36476        1.00     36219     35877
## MSQ4_AG_196 MSQ4_AG_196  49680    49517        1.00     48899     48992
## MSQ4_AG_197 MSQ4_AG_197  47543    47387        1.00     46954     46657
## MSQ4_AG_198 MSQ4_AG_198  30643    30370        0.99     29836     29872
## MSQ4_AG_201 MSQ4_AG_201  19609    19552        1.00     19236     19200
## MSQ4_AG_205 MSQ4_AG_205  27564    27484        1.00     27258     26995
## MSQ4_AG_207 MSQ4_AG_207  47445    47295        1.00     46990     46765
## MSQ4_AG_208 MSQ4_AG_208  49116    48952        1.00     48541     47988
## MSQ4_AG_209 MSQ4_AG_209  47970    47813        1.00     47480     47030
## MSQ4_AG_210 MSQ4_AG_210  45889    45777        1.00     45534     44948
## MSQ4_AG_211 MSQ4_AG_211  33390    33293        1.00     33088     32775
## MSQ4_AG_212 MSQ4_AG_212  42506    42302        1.00     41894     41578
## MSQ4_AG_213 MSQ4_AG_213  43843    43695        1.00     43432     43074
## MSQ4_AG_215 MSQ4_AG_215  52663    52152        0.99     51365     51380
## MSQ4_AG_216 MSQ4_AG_216  52964    52700        1.00     52068     51765
## MSQ4_AG_217 MSQ4_AG_217  47025    46566        0.99     45983     45793
## MSQ4_AG_218 MSQ4_AG_218  52979    52616        0.99     52121     51964
## MSQ4_AG_219 MSQ4_AG_219  41468    41166        0.99     40867     40569
## MSQ4_AG_220 MSQ4_AG_220  47053    46684        0.99     46120     45758
## MSQ4_AG_221 MSQ4_AG_221  48342    48045        0.99     47517     47407
## MSQ4_AG_223 MSQ4_AG_223 177273   176785        1.00    175909    174882
## MSQ4_AG_224 MSQ4_AG_224  60808    60684        1.00     60500     59901
## MSQ4_AG_225 MSQ4_AG_225  64724    64610        1.00     64343     63794
## MSQ4_AG_226 MSQ4_AG_226  65137    64950        1.00     64514     63983
## MSQ4_AG_227 MSQ4_AG_227  50035    49843        1.00     49494     48920
## MSQ4_AG_228 MSQ4_AG_228  47963    47800        1.00     47350     45595
## MSQ4_AG_229 MSQ4_AG_229  51531    51332        1.00     50606     49813
## MSQ4_AG_230 MSQ4_AG_230  47801    47402        0.99     46400     46698
## MSQ4_AG_231 MSQ4_AG_231  47930    47639        0.99     46940     46405
## MSQ4_AG_232 MSQ4_AG_232  53100    52944        1.00     52583     50509
## MSQ4_AG_233 MSQ4_AG_233  48950    48765        1.00     48039     48032
## MSQ4_AG_234 MSQ4_AG_234  53305    53203        1.00     52987     52680
## MSQ4_AG_235 MSQ4_AG_235  48004    47871        1.00     47560     47289
## MSQ4_AG_236 MSQ4_AG_236  46252    46162        1.00     45865     45532
## MSQ4_AG_237 MSQ4_AG_237  44539    44342        1.00     43781     43300
## MSQ4_AG_238 MSQ4_AG_238  36634    36350        0.99     35763     35886
## MSQ4_AG_239 MSQ4_AG_239  53216    52984        1.00     52413     51931
## MSQ4_AG_240 MSQ4_AG_240  56233    56075        1.00     55664     54836
## MSQ4_AG_241 MSQ4_AG_241  51504    51054        0.99     50318     50427
## MSQ4_AG_242 MSQ4_AG_242  60385    60073        0.99     59436     59290
## MSQ4_AG_243 MSQ4_AG_243  37948    37674        0.99     37386     36989
## MSQ4_AG_244 MSQ4_AG_244  54935    54453        0.99     53704     53820
## MSQ4_AG_245 MSQ4_AG_245  46503    46037        0.99     45414     45295
## MSQ4_AG_247 MSQ4_AG_247  54664    54237        0.99     53518     53836
## MSQ4_AG_248 MSQ4_AG_248  55742    55304        0.99     54723     54741
## MSQ4_AG_249 MSQ4_AG_249  49424    49146        0.99     48539     48653
## MSQ6_AG_250 MSQ6_AG_250  44852    44233        0.99     42985     43681
## MSQ6_AG_251 MSQ6_AG_251  72545    72048        0.99     70726     70818
## MSQ6_AG_252 MSQ6_AG_252  76560    76051        0.99     74773     75012
## MSQ6_AG_253 MSQ6_AG_253  62479    62192        1.00     61528     61201
## MSQ6_AG_254 MSQ6_AG_254  78944    78670        1.00     77717     77643
## MSQ6_AG_255 MSQ6_AG_255  72930    72622        1.00     71294     71600
## MSQ6_AG_256 MSQ6_AG_256 194545   193572        0.99    191373    191090
## MSQ6_AG_257 MSQ6_AG_257  62034    61405        0.99     60035     60599
## MSQ6_AG_258 MSQ6_AG_258  77255    76602        0.99     75092     75664
## MSQ6_AG_259 MSQ6_AG_259  58784    58444        0.99     57469     57504
## MSQ6_AG_260 MSQ6_AG_260  65923    65701        1.00     64821     65102
## MSQ6_AG_261 MSQ6_AG_261  51065    50870        1.00     50385     50108
## MSQ6_AG_262 MSQ6_AG_262 123400   123134        1.00    122255    122251
## MSQ6_AG_263 MSQ6_AG_263  60612    60294        0.99     59542     59544
## MSQ6_AG_264 MSQ6_AG_264  66974    66745        1.00     65951     65893
## MSQ6_AG_265 MSQ6_AG_265  59225    58828        0.99     57814     57973
## MSQ6_AG_266 MSQ6_AG_266  51960    51647        0.99     51109     50886
## MSQ6_AG_267 MSQ6_AG_267  46397    46164        0.99     45606     45440
## MSQ6_AG_268 MSQ6_AG_268  58147    57935        1.00     57246     57057
## MSQ6_AG_269 MSQ6_AG_269  39694    39532        1.00     39238     38563
## MSQ6_AG_270 MSQ6_AG_270  68256    68022        1.00     67550     67121
## MSQ6_AG_271 MSQ6_AG_271  58802    58539        1.00     57821     57447
## MSQ6_AG_272 MSQ6_AG_272  60650    60490        1.00     59736     59509
## MSQ6_AG_273 MSQ6_AG_273  57328    57044        1.00     56355     56014
## MSQ6_AG_274 MSQ6_AG_274  29429    29027        0.99     28156     28747
## MSQ6_AG_275 MSQ6_AG_275  66291    65873        0.99     64667     65136
## MSQ6_AG_276 MSQ6_AG_276  72439    71975        0.99     70689     71241
## MSQ6_AG_277 MSQ6_AG_277  54210    53967        1.00     53350     53370
## MSQ6_AG_278 MSQ6_AG_278  53699    53269        0.99     52615     51918
## MSQ6_AG_279 MSQ6_AG_279  57656    57347        0.99     56419     56677
## MSQ6_AG_281 MSQ6_AG_281  74977    74354        0.99     72774     73647
## MSQ6_AG_282 MSQ6_AG_282  46134    45771        0.99     44977     45491
## MSQ6_AG_283 MSQ6_AG_283  43200    42935        0.99     41926     42309
## MSQ6_AG_284 MSQ6_AG_284  61771    61392        0.99     60360     60789
## MSQ6_AG_285 MSQ6_AG_285  52064    51842        1.00     51293     51367
## MSQ6_AG_286 MSQ6_AG_286  60985    60720        1.00     60048     60142
## MSQ6_AG_287 MSQ6_AG_287  58725    58433        1.00     57619     57826
## MSQ6_AG_288 MSQ6_AG_288  53218    52990        1.00     52344     52408
## MSQ6_AG_290 MSQ6_AG_290  35191    34713        0.99     33684     34288
## MSQ6_AG_291 MSQ6_AG_291  59786    59402        0.99     58312     58638
## MSQ6_AG_292 MSQ6_AG_292  71389    70960        0.99     69815     70122
## MSQ6_AG_293 MSQ6_AG_293  48071    47899        1.00     47393     47356
## MSQ6_AG_294 MSQ6_AG_294  60459    60267        1.00     59414     59486
## MSQ6_AG_295 MSQ6_AG_295  56699    56368        0.99     55608     55822
## MSQ6_AG_296 MSQ6_AG_296  52415    52164        1.00     51488     51508
## MSQ6_AG_298 MSQ6_AG_298  40629    40017        0.98     37962     39383
## MSQ6_AG_299 MSQ6_AG_299  67522    67104        0.99     65052     65909
## MSQ6_AG_300 MSQ6_AG_300  62759    62323        0.99     60618     61516
## MSQ6_AG_301 MSQ6_AG_301  51923    51673        1.00     50793     51026
## MSQ6_AG_302 MSQ6_AG_302  67574    67308        1.00     66518     66567
## MSQ6_AG_303 MSQ6_AG_303  69843    69434        0.99     68180     68564
## MSQ6_AG_304 MSQ6_AG_304  53162    52908        1.00     52171     52277
## MSQ6_AG_305 MSQ6_AG_305  89399    88860        0.99     86777     87744
## MSQ6_AG_306 MSQ6_AG_306  47024    46664        0.99     46031     46150
## MSQ6_AG_307 MSQ6_AG_307  52772    52491        0.99     51730     51714
## MSQ6_AG_308 MSQ6_AG_308  57920    57707        1.00     56978     56801
## MSQ6_AG_309 MSQ6_AG_309  45769    45629        1.00     45263     44921
## MSQ6_AG_310 MSQ6_AG_310  69581    69393        1.00     68667     68504
## MSQ6_AG_311 MSQ6_AG_311  51654    51454        1.00     50897     50887
## MSQ6_AG_312 MSQ6_AG_312  47303    47147        1.00     46707     46496
## MSQ6_AG_313 MSQ6_AG_313  59525    59286        1.00     58468     58529
## MSQ6_AG_314 MSQ6_AG_314  57849    57633        1.00     56844     56925
## MSQ6_AG_315 MSQ6_AG_315  66535    66300        1.00     65118     65442
## MSQ6_AG_316 MSQ6_AG_316  55528    55304        1.00     54594     54837
## MSQ6_AG_317 MSQ6_AG_317  44257    44149        1.00     43724     43645
## MSQ6_AG_318 MSQ6_AG_318  68220    68041        1.00     67418     67189
## MSQ6_AG_319 MSQ6_AG_319  67213    66994        1.00     66022     66146
## MSQ6_AG_320 MSQ6_AG_320  48994    48845        1.00     48411     48353
## MSQ6_AG_321 MSQ6_AG_321  64601    64287        1.00     63461     63862
## MSQ6_AG_322 MSQ6_AG_322  44888    44417        0.99     43397     43920
## MSQ6_AG_323 MSQ6_AG_323  52688    52405        0.99     51420     51728
## MSQ6_AG_324 MSQ6_AG_324  51980    51723        1.00     50831     50935
## MSQ6_AG_325 MSQ6_AG_325  47748    47556        1.00     47112     47032
## MSQ6_AG_326 MSQ6_AG_326  56918    56719        1.00     56238     56162
## MSQ6_AG_327 MSQ6_AG_327  47849    47609        0.99     46979     47105
## MSQ6_AG_328 MSQ6_AG_328  50582    50398        1.00     49925     49862
## MSQ6_AG_329 MSQ6_AG_329  45576    45320        0.99     44563     44714
## MSQ6_AG_330 MSQ6_AG_330  37551    37094        0.99     36164     36676
## MSQ6_AG_331 MSQ6_AG_331  49047    48750        0.99     47852     48161
## MSQ6_AG_332 MSQ6_AG_332  55629    55307        0.99     54626     54653
## MSQ6_AG_333 MSQ6_AG_333  39112    38954        1.00     38658     38534
## MSQ6_AG_334 MSQ6_AG_334  56158    55984        1.00     54978     55559
## MSQ6_AG_335 MSQ6_AG_335  56174    55990        1.00     54558     55582
## MSQ6_AG_336 MSQ6_AG_336  49080    48860        1.00     48261     48353
## MSQ6_AG_337 MSQ6_AG_337  51145    50825        0.99     50082     50314
## MSQ6_AG_339 MSQ6_AG_339  61065    60799        1.00     59746     59875
## MSQ6_AG_340 MSQ6_AG_340  49604    49210        0.99     48411     48539
## MSQ6_AG_341 MSQ6_AG_341  39676    39568        1.00     39079     39070
## MSQ6_AG_342 MSQ6_AG_342  67509    67341        1.00     66594     66566
## MSQ6_AG_343 MSQ6_AG_343  54562    54327        1.00     53556     53619
## MSQ6_AG_344 MSQ6_AG_344  51990    51826        1.00     51343     51218
## MSQ6_AG_346 MSQ6_AG_346  56389    56103        0.99     55311     55282
## MSQ6_AG_347 MSQ6_AG_347  52602    52284        0.99     51592     51797
## MSQ6_AG_348 MSQ6_AG_348  55353    55045        0.99     54304     54490
## MSQ6_AG_349 MSQ6_AG_349  47228    47005        1.00     46369     46415
## MSQ6_AG_350 MSQ6_AG_350  51051    50825        1.00     50302     50300
## MSQ6_AG_351 MSQ6_AG_351  50339    50129        1.00     49527     49544
## MSQ6_AG_352 MSQ6_AG_352  53955    53733        1.00     53035     53076
## MSQ6_AG_353 MSQ6_AG_353  50561    50357        1.00     49892     49982
## MSQ6_AG_354 MSQ6_AG_354  55667    55424        1.00     54749     54623
## MSQ6_AG_355 MSQ6_AG_355  53070    52783        0.99     52078     52396
## MSQ6_AG_356 MSQ6_AG_356  43764    43485        0.99     42777     43048
## MSQ6_AG_357 MSQ6_AG_357  55187    55049        1.00     54210     54461
## MSQ6_AG_358 MSQ6_AG_358  44766    44572        1.00     44063     44115
## MSQ6_AG_359 MSQ6_AG_359  60328    60157        1.00     59271     59648
## MSQ6_AG_360 MSQ6_AG_360  58904    58708        1.00     57928     58006
## MSQ6_AG_361 MSQ6_AG_361  45882    45736        1.00     45025     45187
## MSQ6_AG_362 MSQ6_AG_362  53177    53015        1.00     52485     52159
## MSQ6_AG_363 MSQ6_AG_363  37840    37693        1.00     37346     37332
## MSQ6_AG_364 MSQ6_AG_364  48434    48212        1.00     47614     47571
## MSQ6_AG_365 MSQ6_AG_365  42542    42414        1.00     41967     41923
## MSQ6_AG_366 MSQ6_AG_366  50708    50526        1.00     49983     49970
## MSQ6_AG_367 MSQ6_AG_367  47290    47145        1.00     46609     46610
## MSQ6_AG_368 MSQ6_AG_368  53558    53404        1.00     52768     52805
## MSQ6_AG_369 MSQ6_AG_369  48304    48191        1.00     47767     47812
## MSQ6_AG_370 MSQ6_AG_370  50167    49953        1.00     49405     49431
## MSQ6_AG_371 MSQ6_AG_371  48460    48252        1.00     47766     47811
## MSQ6_AG_372 MSQ6_AG_372  51618    51384        1.00     50750     50912
## MSQ6_AG_373 MSQ6_AG_373  38695    38576        1.00     38190     38209
## MSQ6_AG_374 MSQ6_AG_374  31639    31476        0.99     31038     31086
## MSQ6_AG_375 MSQ6_AG_375  57314    57167        1.00     56524     56611
## MSQ6_AG_376 MSQ6_AG_376  26796    26719        1.00     26323     26354
## MSQ6_AG_377 MSQ6_AG_377  42867    42762        1.00     42236     42312
## MSQ6_AG_378 MSQ6_AG_378  49591    49374        1.00     48708     48943
## MSQ6_AG_379 MSQ6_AG_379  49706    49595        1.00     49100     49181
## MSQ6_AG_380 MSQ6_AG_380  82068    81706        1.00     80665     80976
## MSQ6_AG_381 MSQ6_AG_381  52750    52500        1.00     51830     51938
## MSQ6_AG_382 MSQ6_AG_382  73446    73177        1.00     72467     72411
## MSQ6_AG_383 MSQ6_AG_383  69764    69491        1.00     68538     68691
## MSQ6_AG_384 MSQ6_AG_384  61764    61575        1.00     60824     61048
## MSQ6_AG_385 MSQ6_AG_385  52313    52053        1.00     51372     51569
## MSQ6_AG_386 MSQ6_AG_386  51593    51313        0.99     50391     50675
## MSQ6_AG_387 MSQ6_AG_387  46806    46542        0.99     45968     45972
## MSQ6_AG_388 MSQ6_AG_388  47413    47194        1.00     46546     46856
## MSQ6_AG_389 MSQ6_AG_389  58351    58132        1.00     57602     57622
## MSQ6_AG_390 MSQ6_AG_390  49873    49732        1.00     49191     49165
## MSQ6_AG_391 MSQ6_AG_391  63671    63377        1.00     62446     62681
## MSQ6_AG_392 MSQ6_AG_392  49621    49461        1.00     48860     48944
## MSQ6_AG_393 MSQ6_AG_393  47737    47555        1.00     46967     47167
## MSQ6_AG_394 MSQ6_AG_394  61545    61253        1.00     60335     60622
## MSQ6_AG_395 MSQ6_AG_395  53852    53565        0.99     52540     52987
## MSQ6_AG_396 MSQ6_AG_396  64646    64247        0.99     62697     63499
## MSQ6_AG_397 MSQ6_AG_397  48808    48531        0.99     47378     47926
## MSQ6_AG_398 MSQ6_AG_398  47495    47338        1.00     46585     46894
## MSQ6_AG_399 MSQ6_AG_399  87699    87250        0.99     85880     86340
## MSQ6_AG_400 MSQ6_AG_400  70706    70344        0.99     69322     69690
## MSQ6_AG_401 MSQ6_AG_401  60339    60115        1.00     59249     59615
## MSQ6_AG_402 MSQ6_AG_402  50344    50164        1.00     49628     49626
## MSQ6_AG_403 MSQ6_AG_403  44766    44611        1.00     44134     44028
## MSQ6_AG_404 MSQ6_AG_404  45825    45664        1.00     45100     45135
## MSQ6_AG_405 MSQ6_AG_405  49626    49487        1.00     49013     48998
## MSQ6_AG_406 MSQ6_AG_406  57734    57513        1.00     56976     56778
## MSQ6_AG_407 MSQ6_AG_407  85707    85403        1.00     84570     84511
## MSQ6_AG_408 MSQ6_AG_408  72745    72486        1.00     71799     71659
## MSQ6_AG_409 MSQ6_AG_409  50870    50700        1.00     50186     50147
## MSQ6_AG_410 MSQ6_AG_410  44125    43923        1.00     43307     43366
## MSQ6_AG_411 MSQ6_AG_411  37214    37112        1.00     36637     36662
## MSQ6_AG_412 MSQ6_AG_412  57349    57150        1.00     56372     56668
## MSQ6_AG_413 MSQ6_AG_413  44970    44835        1.00     44400     44409
## MSQ6_AG_414 MSQ6_AG_414  43763    43641        1.00     43149     43180
## MSQ6_AG_415 MSQ6_AG_415  60177    60012        1.00     59247     59381
## MSQ6_AG_416 MSQ6_AG_416  51357    51261        1.00     50857     50884
## MSQ6_AG_417 MSQ6_AG_417  41218    41094        1.00     40645     40533
## MSQ6_AG_418 MSQ6_AG_418  53032    52842        1.00     52108     52151
## MSQ6_AG_419 MSQ6_AG_419  45673    45546        1.00     45060     45092
## MSQ6_AG_420 MSQ6_AG_420  46538    46392        1.00     45733     45945
## MSQ6_AG_421 MSQ6_AG_421  47359    47215        1.00     46592     46627
## MSQ6_AG_422 MSQ6_AG_422  56935    56791        1.00     56147     56148
## MSQ6_AG_423 MSQ6_AG_423  51078    50931        1.00     50170     50260
## MSQ6_AG_424 MSQ6_AG_424  43618    43500        1.00     42926     42941
## MSQ6_AG_425 MSQ6_AG_425  47220    47071        1.00     46403     46542
## MSQ6_AG_426 MSQ6_AG_426  54315    54183        1.00     53607     53647
## MSQ6_AG_427 MSQ6_AG_427  54083    53914        1.00     53211     53282
## MSQ6_AG_428 MSQ6_AG_428  59050    58793        1.00     57902     58192
## MSQ6_AG_429 MSQ6_AG_429  58474    58254        1.00     57537     57642
## MSQ6_AG_430 MSQ6_AG_430  51572    51402        1.00     50774     50701
## MSQ6_AG_431 MSQ6_AG_431  59620    59481        1.00     58897     59038
## MSQ6_AG_432 MSQ6_AG_432  60644    60467        1.00     59646     59775
## MSQ6_AG_433 MSQ6_AG_433  54252    54048        1.00     53232     53335
## MSQ6_AG_434 MSQ6_AG_434  47447    47237        1.00     46574     46630
## MSQ6_AG_490 MSQ6_AG_490  66744    66455        1.00     65323     65479
## MSQ6_AG_491 MSQ6_AG_491  65688    65409        1.00     64553     64591
## MSQ6_AG_492 MSQ6_AG_492  42282    42132        1.00     41655     41507
## MSQ6_AG_493 MSQ6_AG_493  70679    70422        1.00     69786     69632
## MSQ6_AG_494 MSQ6_AG_494  64149    63878        1.00     62852     63031
## MSQ6_AG_495 MSQ6_AG_495  66106    65772        0.99     65139     65014
## MSQ6_AG_496 MSQ6_AG_496  67344    67025        1.00     65804     66358
## MSQ6_AG_497 MSQ6_AG_497  47094    46707        0.99     45473     46244
## MSQ6_AG_498 MSQ6_AG_498  67572    67267        1.00     65981     66428
## MSQ6_AG_499 MSQ6_AG_499  67255    66994        1.00     65937     66416
## MSQ6_AG_500 MSQ6_AG_500  57626    57445        1.00     56882     56862
## MSQ6_AG_501 MSQ6_AG_501  69800    69629        1.00     69042     68944
## MSQ6_AG_502 MSQ6_AG_502  66731    66475        1.00     65494     65868
## MSQ6_AG_503 MSQ6_AG_503  45499    45355        1.00     44891     45065
## MSQ6_AG_504 MSQ6_AG_504  79099    78602        0.99     77137     77843
## MSQ6_AG_505 MSQ6_AG_505  19457    19274        0.99     18571     18901
## MSQ6_AG_506 MSQ6_AG_506  24362    24239        0.99     23778     23869
## MSQ6_AG_509 MSQ6_AG_509  19742    19663        1.00     19494     19373
## MSQ6_AG_511 MSQ6_AG_511  27149    27010        0.99     26817     26739
## MSQ6_AG_512 MSQ6_AG_512  29768    29635        1.00     29111     29208
## MSQ6_AG_514 MSQ6_AG_514  76385    75928        0.99     74202     74970
## MSQ6_AG_515 MSQ6_AG_515  80169    79812        1.00     78438     79046
## MSQ6_AG_516 MSQ6_AG_516  59685    59446        1.00     58799     58725
## MSQ6_AG_517 MSQ6_AG_517  69793    69596        1.00     68933     69150
## MSQ6_AG_518 MSQ6_AG_518  73792    73415        0.99     72195     72677
## MSQ6_AG_519 MSQ6_AG_519  64707    64401        1.00     63506     63635
## MSQ6_AG_520 MSQ6_AG_520  69869    69411        0.99     67836     68545
## MSQ6_AG_521 MSQ6_AG_521  22876    22575        0.99     21552     22319
## MSQ6_AG_522 MSQ6_AG_522  70853    70436        0.99     69010     69691
## MSQ6_AG_523 MSQ6_AG_523  66286    66030        1.00     64938     65385
## MSQ6_AG_524 MSQ6_AG_524  63794    63539        1.00     62844     62818
## MSQ6_AG_525 MSQ6_AG_525  76951    76627        1.00     76016     76028
## MSQ6_AG_526 MSQ6_AG_526  61480    61249        1.00     60023     60475
## MSQ6_AG_527 MSQ6_AG_527  65287    65039        1.00     64261     64598
## MSQ6_AG_528 MSQ6_AG_528  62397    61982        0.99     60877     61304
## MSQ6_AG_529 MSQ6_AG_529  40902    40453        0.99     39440     40103
## MSQ6_AG_530 MSQ6_AG_530  67854    67566        1.00     66284     66723
##             denoisedF_pc denoisedR_pc merged merged_pc tabled chimera_out
## MSQ4_AG_012         0.99         0.99  84118      0.91  84118       63132
## MSQ4_AG_013         0.99         0.99  35885      0.89  35885       28004
## MSQ4_AG_014         0.99         0.99  47579      0.91  47579       34919
## MSQ4_AG_015         0.99         0.98  50581      0.92  50581       36310
## MSQ4_AG_016         0.98         0.98  34271      0.87  34271       25798
## MSQ4_AG_017         0.99         0.99  44660      0.87  44660       31556
## MSQ4_AG_018         0.99         0.99  45468      0.92  45468       39875
## MSQ4_AG_019         0.99         0.99  44110      0.88  44110       33557
## MSQ4_AG_020         0.99         0.99  62578      0.91  62578       43716
## MSQ4_AG_021         0.99         0.99  29195      0.84  29195       23176
## MSQ4_AG_022         0.99         0.99  50063      0.91  50063       41438
## MSQ4_AG_023         0.99         0.99  48417      0.92  48417       40845
## MSQ4_AG_024         0.99         0.99  40356      0.88  40356       28168
## MSQ4_AG_025         0.99         0.99  51495      0.90  51495       38257
## MSQ4_AG_027         0.99         0.99  53208      0.90  53208       36624
## MSQ4_AG_028         0.99         0.98  37204      0.85  37204       26811
## MSQ4_AG_029         0.99         0.99  34583      0.89  34583       24267
## MSQ4_AG_030         0.99         0.99  62905      0.85  62905       47010
## MSQ4_AG_031         0.98         0.98  48893      0.87  48893       38050
## MSQ4_AG_033         0.99         0.99  64489      0.91  64489       44462
## MSQ4_AG_034         0.99         0.99  46671      0.89  46671       33723
## MSQ4_AG_035         0.99         0.99  59596      0.92  59596       48208
## MSQ4_AG_036         0.99         0.98  67973      0.89  67973       47982
## MSQ4_AG_037         0.99         0.99  29945      0.82  29945       21250
## MSQ4_AG_038         0.99         0.99  46475      0.89  46475       34480
## MSQ4_AG_039         0.98         0.99  40530      0.85  40530       28166
## MSQ4_AG_040         0.98         0.99  31959      0.85  31959       23637
## MSQ4_AG_041         0.98         0.99  39614      0.90  39614       28883
## MSQ4_AG_042         0.99         0.99  47012      0.90  47012       36332
## MSQ4_AG_043         0.99         0.99  46794      0.92  46794       36253
## MSQ4_AG_044         0.99         0.99  45721      0.88  45721       33676
## MSQ4_AG_045         0.99         0.99  33004      0.86  33004       22118
## MSQ4_AG_046         0.98         0.99  41812      0.84  41812       29882
## MSQ4_AG_047         0.98         0.98  45976      0.92  45976       35244
## MSQ4_AG_049         0.98         0.98  52957      0.91  52957       37365
## MSQ4_AG_050         0.98         0.98  44788      0.87  44788       33113
## MSQ4_AG_051         0.98         0.99  48832      0.86  48832       33451
## MSQ4_AG_052         0.99         0.99  56975      0.91  56975       43344
## MSQ4_AG_053         0.98         0.99  38222      0.86  38222       26692
## MSQ4_AG_054         0.98         0.99  53631      0.88  53631       35882
## MSQ4_AG_055         0.98         0.99  50266      0.83  50266       34369
## MSQ4_AG_056         0.95         0.98  27651      0.85  27651       18872
## MSQ4_AG_057         0.99         0.99  59495      0.95  59495       53231
## MSQ4_AG_058         0.99         0.98  55812      0.94  55812       46387
## MSQ4_AG_059         0.99         0.99  57344      0.95  57344       52456
## MSQ4_AG_062         0.99         0.98  85117      0.95  85117       78943
## MSQ4_AG_064         0.98         0.98  44827      0.94  44827       42321
## MSQ4_AG_065         0.98         0.98  47461      0.96  47461       43474
## MSQ4_AG_072         0.99         0.98  26150      0.95  26150       24746
## MSQ4_AG_076         0.99         0.98  58443      0.95  58443       55585
## MSQ4_AG_077         0.99         0.95  42985      0.91  42985       37756
## MSQ4_AG_080         0.98         0.98  38686      0.94  38686       36761
## MSQ4_AG_083         0.99         0.98  50942      0.93  50942       39045
## MSQ4_AG_084         0.99         0.99  59641      0.94  59641       45763
## MSQ4_AG_085         0.99         0.98  42947      0.90  42947       31589
## MSQ4_AG_086         0.99         0.98  61486      0.90  61486       47584
## MSQ4_AG_087         0.98         0.99  46439      0.95  46439       40651
## MSQ4_AG_088         0.99         0.99  36949      0.94  36949       32150
## MSQ4_AG_089         0.99         0.99  54586      0.94  54586       45583
## MSQ4_AG_091         0.99         0.99  57767      0.98  57767       57555
## MSQ4_AG_092         0.99         0.99  63673      0.96  63673       50090
## MSQ4_AG_093         0.99         0.99  46470      0.94  46470       35237
## MSQ4_AG_095         0.99         0.98  47906      0.91  47906       37534
## MSQ4_AG_096         0.99         0.99  53595      0.93  53595       34866
## MSQ4_AG_101         0.99         0.98  39584      0.90  39584       32686
## MSQ4_AG_103         0.99         0.98  38895      0.91  38895       30485
## MSQ4_AG_104         0.99         0.98  38210      0.92  38210       28858
## MSQ4_AG_105         0.99         0.98  40867      0.92  40867       30798
## MSQ4_AG_106         0.99         0.98  24486      0.90  24486       18531
## MSQ4_AG_107         0.99         0.98  37138      0.89  37138       26957
## MSQ4_AG_108         0.99         0.98  35180      0.91  35180       27091
## MSQ4_AG_109         0.99         0.98  37843      0.91  37843       29029
## MSQ4_AG_111         0.99         0.97  36671      0.91  36671       26044
## MSQ4_AG_112         0.99         0.99  50090      0.95  50090       45944
## MSQ4_AG_113         0.99         0.98  43348      0.92  43348       36643
## MSQ4_AG_117         0.99         0.98  41168      0.89  41168       30639
## MSQ4_AG_119         0.99         0.97  35402      0.89  35402       27533
## MSQ4_AG_121         0.99         0.98  49842      0.91  49842       35845
## MSQ4_AG_125         0.99         0.99  36344      0.94  36344       29264
## MSQ4_AG_127         0.99         0.98  35583      0.89  35583       27454
## MSQ4_AG_128         0.99         0.98  57515      0.93  57515       45962
## MSQ4_AG_129         0.99         0.98  47481      0.93  47481       34436
## MSQ4_AG_130         0.99         0.99  39997      0.91  39997       30300
## MSQ4_AG_131         0.99         0.98  42324      0.92  42324       28664
## MSQ4_AG_132         0.99         0.99  56666      0.92  56666       44270
## MSQ4_AG_133         1.00         1.00  51772      0.98  51772       47432
## MSQ4_AG_135         0.99         0.99  44957      0.95  44957       39379
## MSQ4_AG_136         0.99         0.99  45705      0.92  45705       34867
## MSQ4_AG_141         0.99         0.98  36976      0.88  36976       28476
## MSQ4_AG_143         0.99         0.99  46303      0.90  46303       35175
## MSQ4_AG_144         0.99         0.98  51711      0.93  51711       35801
## MSQ4_AG_145         0.99         0.97  54728      0.90  54728       36833
## MSQ4_AG_148         0.98         0.98  94505      0.91  94505       68895
## MSQ4_AG_149         0.99         0.99  56970      0.88  56970       41931
## MSQ4_AG_151         0.98         0.99  43113      0.85  43113       31350
## MSQ4_AG_152         0.99         0.98  53292      0.90  53292       39031
## MSQ4_AG_153         0.99         0.99  66121      0.91  66121       49753
## MSQ4_AG_154         0.99         0.99  51921      0.88  51921       38218
## MSQ4_AG_155         0.99         0.98  45499      0.90  45499       30522
## MSQ4_AG_156         0.99         0.99  44992      0.89  44992       33495
## MSQ4_AG_157         0.99         0.99  45726      0.88  45726       34703
## MSQ4_AG_159         0.99         0.99  44664      0.90  44664       34623
## MSQ4_AG_160         0.99         0.99  46331      0.90  46331       34466
## MSQ4_AG_161         0.99         0.99  47183      0.90  47183       36003
## MSQ4_AG_162         0.99         0.99  24387      0.92  24387       18381
## MSQ4_AG_163         0.99         0.99  35645      0.88  35645       26906
## MSQ4_AG_164         0.99         0.99  52609      0.91  52609       40076
## MSQ4_AG_165         0.99         0.99  40818      0.94  40818       31663
## MSQ4_AG_167         0.99         0.98  36059      0.90  36059       27235
## MSQ4_AG_168         0.99         0.99  53278      0.90  53278       38492
## MSQ4_AG_169         0.99         0.99  46600      0.90  46600       35599
## MSQ4_AG_170         0.99         0.99  46960      0.90  46960       34791
## MSQ4_AG_171         0.99         0.99  49114      0.90  49114       34171
## MSQ4_AG_172         0.99         0.99  44662      0.90  44662       33548
## MSQ4_AG_173         0.99         0.99  47078      0.90  47078       33771
## MSQ4_AG_175         0.99         0.99  56264      0.89  56264       37348
## MSQ4_AG_176         0.99         0.99  48469      0.93  48469       38771
## MSQ4_AG_177         0.99         0.99  39913      0.90  39913       29662
## MSQ4_AG_178         0.99         0.99  35333      0.91  35333       24028
## MSQ4_AG_179         0.99         0.99  36750      0.91  36750       27466
## MSQ4_AG_180         0.99         0.99  39455      0.95  39455       31624
## MSQ4_AG_181         0.99         0.99  36247      0.89  36247       28373
## MSQ4_AG_183         0.99         0.99  36842      0.91  36842       28511
## MSQ4_AG_184         0.99         0.99  42946      0.88  42946       30442
## MSQ4_AG_185         0.99         0.99  42227      0.88  42227       31339
## MSQ4_AG_186         0.99         0.99  42245      0.90  42245       31727
## MSQ4_AG_187         0.99         0.98  31145      0.85  31145       23085
## MSQ4_AG_188         0.99         0.99  34331      0.88  34331       25966
## MSQ4_AG_189         0.99         0.99  48798      0.89  48798       34955
## MSQ4_AG_191         0.98         0.99  51412      0.89  51412       38217
## MSQ4_AG_192         0.99         0.99  43926      0.95  43926       37249
## MSQ4_AG_193         0.99         0.99  44846      0.95  44846       39378
## MSQ4_AG_194         1.00         0.99  41202      0.93  41202       32588
## MSQ4_AG_195         0.99         0.98  33188      0.92  33188       24980
## MSQ4_AG_196         0.99         0.99  44186      0.90  44186       32211
## MSQ4_AG_197         0.99         0.98  42794      0.91  42794       31808
## MSQ4_AG_198         0.98         0.98  27589      0.92  27589       21792
## MSQ4_AG_201         0.98         0.98  17017      0.88  17017       12911
## MSQ4_AG_205         0.99         0.98  25057      0.92  25057       22030
## MSQ4_AG_207         0.99         0.99  44325      0.94  44325       35316
## MSQ4_AG_208         0.99         0.98  43661      0.90  43661       31211
## MSQ4_AG_209         0.99         0.98  42926      0.90  42926       31033
## MSQ4_AG_210         0.99         0.98  41208      0.90  41208       32410
## MSQ4_AG_211         0.99         0.98  30328      0.92  30328       22976
## MSQ4_AG_212         0.99         0.98  37441      0.89  37441       28577
## MSQ4_AG_213         0.99         0.99  40164      0.92  40164       30989
## MSQ4_AG_215         0.98         0.99  45323      0.88  45323       29770
## MSQ4_AG_216         0.99         0.98  45890      0.88  45890       31999
## MSQ4_AG_217         0.99         0.98  41446      0.90  41446       25845
## MSQ4_AG_218         0.99         0.99  46997      0.90  46997       29813
## MSQ4_AG_219         0.99         0.99  37261      0.91  37261       29063
## MSQ4_AG_220         0.99         0.98  41383      0.90  41383       27283
## MSQ4_AG_221         0.99         0.99  42561      0.90  42561       28554
## MSQ4_AG_223         1.00         0.99 165714      0.94 165714      121734
## MSQ4_AG_224         1.00         0.99  56715      0.94  56715       43308
## MSQ4_AG_225         1.00         0.99  60923      0.95  60923       44919
## MSQ4_AG_226         0.99         0.99  58760      0.91  58760       42654
## MSQ4_AG_227         0.99         0.98  44646      0.90  44646       30151
## MSQ4_AG_228         0.99         0.95  35223      0.74  35223       26182
## MSQ4_AG_229         0.99         0.97  43573      0.86  43573       34084
## MSQ4_AG_230         0.98         0.99  41460      0.89  41460       27509
## MSQ4_AG_231         0.99         0.97  40212      0.86  40212       29257
## MSQ4_AG_232         0.99         0.95  45053      0.86  45053       39903
## MSQ4_AG_233         0.99         0.98  43450      0.90  43450       31217
## MSQ4_AG_234         1.00         0.99  50677      0.96  50677       36152
## MSQ4_AG_235         0.99         0.99  44432      0.93  44432       34859
## MSQ4_AG_236         0.99         0.99  43005      0.94  43005       35312
## MSQ4_AG_237         0.99         0.98  37957      0.87  37957       27057
## MSQ4_AG_238         0.98         0.99  32481      0.91  32481       25671
## MSQ4_AG_239         0.99         0.98  47954      0.91  47954       38863
## MSQ4_AG_240         0.99         0.98  49795      0.89  49795       37376
## MSQ4_AG_241         0.99         0.99  44402      0.88  44402       28959
## MSQ4_AG_242         0.99         0.99  53955      0.91  53955       36918
## MSQ4_AG_243         0.99         0.98  33123      0.89  33123       22034
## MSQ4_AG_244         0.99         0.99  48755      0.91  48755       35309
## MSQ4_AG_245         0.99         0.98  40621      0.89  40621       25726
## MSQ4_AG_247         0.99         0.99  49391      0.92  49391       39348
## MSQ4_AG_248         0.99         0.99  50107      0.92  50107       33598
## MSQ4_AG_249         0.99         0.99  44587      0.92  44587       33624
## MSQ6_AG_250         0.97         0.99  37694      0.88  37694       26929
## MSQ6_AG_251         0.98         0.98  63881      0.90  63881       44125
## MSQ6_AG_252         0.98         0.99  67904      0.91  67904       49814
## MSQ6_AG_253         0.99         0.98  54905      0.89  54905       42532
## MSQ6_AG_254         0.99         0.99  68254      0.88  68254       49859
## MSQ6_AG_255         0.98         0.99  60502      0.85  60502       44187
## MSQ6_AG_256         0.99         0.99 173270      0.91 173270      124781
## MSQ6_AG_257         0.98         0.99  53780      0.90  53780       39110
## MSQ6_AG_258         0.98         0.99  65965      0.88  65965       47209
## MSQ6_AG_259         0.98         0.98  52534      0.91  52534       39238
## MSQ6_AG_260         0.99         0.99  59729      0.92  59729       49569
## MSQ6_AG_261         0.99         0.99  45226      0.90  45226       32999
## MSQ6_AG_262         0.99         0.99 113946      0.93 113946       89836
## MSQ6_AG_263         0.99         0.99  53809      0.90  53809       37806
## MSQ6_AG_264         0.99         0.99  58534      0.89  58534       44726
## MSQ6_AG_265         0.98         0.99  52635      0.91  52635       37761
## MSQ6_AG_266         0.99         0.99  45695      0.89  45695       32968
## MSQ6_AG_267         0.99         0.98  40905      0.90  40905       29158
## MSQ6_AG_268         0.99         0.98  50818      0.89  50818       37290
## MSQ6_AG_269         0.99         0.98  33943      0.87  33943       24804
## MSQ6_AG_270         0.99         0.99  61876      0.92  61876       45108
## MSQ6_AG_271         0.99         0.98  52382      0.91  52382       35201
## MSQ6_AG_272         0.99         0.98  51500      0.86  51500       37456
## MSQ6_AG_273         0.99         0.98  50104      0.89  50104       33862
## MSQ6_AG_274         0.97         0.99  24422      0.87  24422       18014
## MSQ6_AG_275         0.98         0.99  59830      0.93  59830       43190
## MSQ6_AG_276         0.98         0.99  65709      0.93  65709       45910
## MSQ6_AG_277         0.99         0.99  49868      0.93  49868       36000
## MSQ6_AG_278         0.99         0.97  47465      0.90  47465       34478
## MSQ6_AG_279         0.98         0.99  51758      0.92  51758       36188
## MSQ6_AG_281         0.98         0.99  68075      0.94  68075       50265
## MSQ6_AG_282         0.98         0.99  42248      0.94  42248       36121
## MSQ6_AG_283         0.98         0.99  37345      0.89  37345       27664
## MSQ6_AG_284         0.98         0.99  55612      0.92  55612       39123
## MSQ6_AG_285         0.99         0.99  47517      0.93  47517       34414
## MSQ6_AG_286         0.99         0.99  55956      0.93  55956       41213
## MSQ6_AG_287         0.99         0.99  53508      0.93  53508       39884
## MSQ6_AG_288         0.99         0.99  48415      0.92  48415       35603
## MSQ6_AG_290         0.97         0.99  29649      0.88  29649       22393
## MSQ6_AG_291         0.98         0.99  53165      0.91  53165       37401
## MSQ6_AG_292         0.98         0.99  63752      0.91  63752       44075
## MSQ6_AG_293         0.99         0.99  42784      0.90  42784       33102
## MSQ6_AG_294         0.99         0.99  53184      0.90  53184       39177
## MSQ6_AG_295         0.99         0.99  51337      0.92  51337       36374
## MSQ6_AG_296         0.99         0.99  47000      0.91  47000       33074
## MSQ6_AG_298         0.95         0.98  33276      0.88  33276       25918
## MSQ6_AG_299         0.97         0.98  58828      0.90  58828       44286
## MSQ6_AG_300         0.97         0.99  55474      0.92  55474       41040
## MSQ6_AG_301         0.98         0.99  47109      0.93  47109       33787
## MSQ6_AG_302         0.99         0.99  62749      0.94  62749       49311
## MSQ6_AG_303         0.98         0.99  63342      0.93  63342       47625
## MSQ6_AG_304         0.99         0.99  48976      0.94  48976       38035
## MSQ6_AG_305         0.98         0.99  80614      0.93  80614       59937
## MSQ6_AG_306         0.99         0.99  42485      0.92  42485       30890
## MSQ6_AG_307         0.99         0.99  47698      0.92  47698       36986
## MSQ6_AG_308         0.99         0.98  52627      0.92  52627       40691
## MSQ6_AG_309         0.99         0.98  41374      0.91  41374       31897
## MSQ6_AG_310         0.99         0.99  64030      0.93  64030       47728
## MSQ6_AG_311         0.99         0.99  47420      0.93  47420       37306
## MSQ6_AG_312         0.99         0.99  43385      0.93  43385       32803
## MSQ6_AG_313         0.99         0.99  53678      0.92  53678       39652
## MSQ6_AG_314         0.99         0.99  52397      0.92  52397       39270
## MSQ6_AG_315         0.98         0.99  60172      0.92  60172       46336
## MSQ6_AG_316         0.99         0.99  52252      0.96  52252       45015
## MSQ6_AG_317         0.99         0.99  40660      0.93  40660       30947
## MSQ6_AG_318         0.99         0.99  62738      0.93  62738       46795
## MSQ6_AG_319         0.99         0.99  60641      0.92  60641       43420
## MSQ6_AG_320         0.99         0.99  45847      0.95  45847       36090
## MSQ6_AG_321         0.99         0.99  61290      0.97  61290       49278
## MSQ6_AG_322         0.98         0.99  39771      0.92  39771       30670
## MSQ6_AG_323         0.98         0.99  47306      0.92  47306       36889
## MSQ6_AG_324         0.98         0.98  48682      0.96  48682       40654
## MSQ6_AG_325         0.99         0.99  44695      0.95  44695       35242
## MSQ6_AG_326         0.99         0.99  53216      0.95  53216       44122
## MSQ6_AG_327         0.99         0.99  44790      0.95  44790       35011
## MSQ6_AG_328         0.99         0.99  47799      0.96  47799       36199
## MSQ6_AG_329         0.98         0.99  41344      0.93  41344       32872
## MSQ6_AG_330         0.97         0.99  32859      0.91  32859       25123
## MSQ6_AG_331         0.98         0.99  44429      0.93  44429       34215
## MSQ6_AG_332         0.99         0.99  50825      0.93  50825       39342
## MSQ6_AG_333         0.99         0.99  36206      0.94  36206       27584
## MSQ6_AG_334         0.98         0.99  51798      0.94  51798       36753
## MSQ6_AG_335         0.97         0.99  50705      0.93  50705       34090
## MSQ6_AG_336         0.99         0.99  43798      0.91  43798       32369
## MSQ6_AG_337         0.99         0.99  46319      0.92  46319       33486
## MSQ6_AG_339         0.98         0.98  51026      0.85  51026       37026
## MSQ6_AG_340         0.98         0.99  43102      0.89  43102       30284
## MSQ6_AG_341         0.99         0.99  33431      0.86  33431       24770
## MSQ6_AG_342         0.99         0.99  58799      0.88  58799       41706
## MSQ6_AG_343         0.99         0.99  47333      0.88  47333       34794
## MSQ6_AG_344         0.99         0.99  48093      0.94  48093       38160
## MSQ6_AG_346         0.99         0.99  50234      0.91  50234       37273
## MSQ6_AG_347         0.99         0.99  48036      0.93  48036       35747
## MSQ6_AG_348         0.99         0.99  50207      0.92  50207       36364
## MSQ6_AG_349         0.99         0.99  43142      0.93  43142       32689
## MSQ6_AG_350         0.99         0.99  46711      0.93  46711       35016
## MSQ6_AG_351         0.99         0.99  45769      0.92  45769       33891
## MSQ6_AG_352         0.99         0.99  48984      0.92  48984       35926
## MSQ6_AG_353         0.99         0.99  47055      0.94  47055       33868
## MSQ6_AG_354         0.99         0.99  50335      0.92  50335       38224
## MSQ6_AG_355         0.99         0.99  49621      0.95  49621       41368
## MSQ6_AG_356         0.98         0.99  39705      0.93  39705       29213
## MSQ6_AG_357         0.98         0.99  49886      0.92  49886       37943
## MSQ6_AG_358         0.99         0.99  41415      0.94  41415       33353
## MSQ6_AG_359         0.99         0.99  55103      0.93  55103       41571
## MSQ6_AG_360         0.99         0.99  53360      0.92  53360       39599
## MSQ6_AG_361         0.98         0.99  41682      0.93  41682       32148
## MSQ6_AG_362         0.99         0.98  46406      0.88  46406       35602
## MSQ6_AG_363         0.99         0.99  35087      0.94  35087       27164
## MSQ6_AG_364         0.99         0.99  43897      0.92  43897       35941
## MSQ6_AG_365         0.99         0.99  39240      0.94  39240       30861
## MSQ6_AG_366         0.99         0.99  46612      0.93  46612       37278
## MSQ6_AG_367         0.99         0.99  43131      0.93  43131       34694
## MSQ6_AG_368         0.99         0.99  49652      0.94  49652       37761
## MSQ6_AG_369         0.99         0.99  46123      0.97  46123       38432
## MSQ6_AG_370         0.99         0.99  46317      0.94  46317       35563
## MSQ6_AG_371         0.99         0.99  45118      0.94  45118       35425
## MSQ6_AG_372         0.99         0.99  47667      0.94  47667       37388
## MSQ6_AG_373         0.99         0.99  36114      0.95  36114       30253
## MSQ6_AG_374         0.99         0.99  28431      0.92  28431       22578
## MSQ6_AG_375         0.99         0.99  52379      0.93  52379       41121
## MSQ6_AG_376         0.99         0.99  23818      0.90  23818       18120
## MSQ6_AG_377         0.99         0.99  38771      0.92  38771       30174
## MSQ6_AG_378         0.99         0.99  45983      0.94  45983       34059
## MSQ6_AG_379         0.99         0.99  46988      0.96  46988       36536
## MSQ6_AG_380         0.99         0.99  75663      0.94  75663       54422
## MSQ6_AG_381         0.99         0.99  48047      0.93  48047       33745
## MSQ6_AG_382         0.99         0.99  67645      0.93  67645       47662
## MSQ6_AG_383         0.99         0.99  63239      0.92  63239       45958
## MSQ6_AG_384         0.99         0.99  56477      0.93  56477       45184
## MSQ6_AG_385         0.99         0.99  48051      0.94  48051       34946
## MSQ6_AG_386         0.98         0.99  46169      0.92  46169       32861
## MSQ6_AG_387         0.99         0.99  41514      0.90  41514       29841
## MSQ6_AG_388         0.99         0.99  43339      0.93  43339       32762
## MSQ6_AG_389         0.99         0.99  54221      0.94  54221       42292
## MSQ6_AG_390         0.99         0.99  45346      0.92  45346       34330
## MSQ6_AG_391         0.99         0.99  58233      0.93  58233       45537
## MSQ6_AG_392         0.99         0.99  45370      0.93  45370       34869
## MSQ6_AG_393         0.99         0.99  44005      0.94  44005       34710
## MSQ6_AG_394         0.99         0.99  55545      0.92  55545       42780
## MSQ6_AG_395         0.98         0.99  48288      0.92  48288       35991
## MSQ6_AG_396         0.98         0.99  58092      0.93  58092       45382
## MSQ6_AG_397         0.98         0.99  43474      0.92  43474       34890
## MSQ6_AG_398         0.98         0.99  42330      0.91  42330       32099
## MSQ6_AG_399         0.98         0.99  79718      0.93  79718       60535
## MSQ6_AG_400         0.99         0.99  64651      0.93  64651       49432
## MSQ6_AG_401         0.99         0.99  55019      0.93  55019       41237
## MSQ6_AG_402         0.99         0.99  46033      0.93  46033       35908
## MSQ6_AG_403         0.99         0.99  39557      0.90  39557       29153
## MSQ6_AG_404         0.99         0.99  41748      0.93  41748       31418
## MSQ6_AG_405         0.99         0.99  45464      0.93  45464       35131
## MSQ6_AG_406         0.99         0.99  52407      0.92  52407       41200
## MSQ6_AG_407         0.99         0.99  78971      0.93  78971       61652
## MSQ6_AG_408         0.99         0.99  66877      0.93  66877       51398
## MSQ6_AG_409         0.99         0.99  46648      0.93  46648       36318
## MSQ6_AG_410         0.99         0.99  40215      0.93  40215       31866
## MSQ6_AG_411         0.99         0.99  33332      0.91  33332       24987
## MSQ6_AG_412         0.99         0.99  52580      0.93  52580       40913
## MSQ6_AG_413         0.99         0.99  41307      0.93  41307       33713
## MSQ6_AG_414         0.99         0.99  40248      0.93  40248       31111
## MSQ6_AG_415         0.99         0.99  55199      0.93  55199       42346
## MSQ6_AG_416         0.99         0.99  48857      0.96  48857       37759
## MSQ6_AG_417         0.99         0.99  37492      0.92  37492       28740
## MSQ6_AG_418         0.99         0.99  47285      0.91  47285       36412
## MSQ6_AG_419         0.99         0.99  40991      0.91  40991       30021
## MSQ6_AG_420         0.99         0.99  42065      0.92  42065       32891
## MSQ6_AG_421         0.99         0.99  41828      0.90  41828       32405
## MSQ6_AG_422         0.99         0.99  51354      0.91  51354       40869
## MSQ6_AG_423         0.99         0.99  45144      0.90  45144       33271
## MSQ6_AG_424         0.99         0.99  38971      0.91  38971       29615
## MSQ6_AG_425         0.99         0.99  42221      0.91  42221       31675
## MSQ6_AG_426         0.99         0.99  50675      0.95  50675       37581
## MSQ6_AG_427         0.99         0.99  48190      0.91  48190       36710
## MSQ6_AG_428         0.98         0.99  52379      0.90  52379       39430
## MSQ6_AG_429         0.99         0.99  51657      0.90  51657       40730
## MSQ6_AG_430         0.99         0.99  44920      0.88  44920       34324
## MSQ6_AG_431         0.99         0.99  55504      0.94  55504       47958
## MSQ6_AG_432         0.99         0.99  53588      0.90  53588       42098
## MSQ6_AG_433         0.98         0.99  47239      0.89  47239       36215
## MSQ6_AG_434         0.99         0.99  40830      0.88  40830       30546
## MSQ6_AG_490         0.98         0.99  60561      0.93  60561       52584
## MSQ6_AG_491         0.99         0.99  60150      0.93  60150       50751
## MSQ6_AG_492         0.99         0.99  38560      0.93  38560       33283
## MSQ6_AG_493         0.99         0.99  66477      0.95  66477       58699
## MSQ6_AG_494         0.98         0.99  60377      0.96  60377       58881
## MSQ6_AG_495         0.99         0.99  61479      0.94  61479       50758
## MSQ6_AG_496         0.98         0.99  63124      0.96  63124       60863
## MSQ6_AG_497         0.97         0.99  42723      0.94  42723       36375
## MSQ6_AG_498         0.98         0.99  61694      0.94  61694       50918
## MSQ6_AG_499         0.98         0.99  63155      0.96  63155       58042
## MSQ6_AG_500         0.99         0.99  54043      0.95  54043       46975
## MSQ6_AG_501         0.99         0.99  65940      0.96  65940       57531
## MSQ6_AG_502         0.99         0.99  62165      0.95  62165       53894
## MSQ6_AG_503         0.99         0.99  44053      0.98  44053       43444
## MSQ6_AG_504         0.98         0.99  73581      0.95  73581       63433
## MSQ6_AG_505         0.96         0.98  17045      0.92  17045       14674
## MSQ6_AG_506         0.98         0.98  22226      0.93  22226       19782
## MSQ6_AG_509         0.99         0.99  18377      0.94  18377       15767
## MSQ6_AG_511         0.99         0.99  25587      0.95  25587       22168
## MSQ6_AG_512         0.98         0.99  26676      0.92  26676       22917
## MSQ6_AG_514         0.98         0.99  69097      0.93  69097       56620
## MSQ6_AG_515         0.98         0.99  75085      0.96  75085       68258
## MSQ6_AG_516         0.99         0.99  55577      0.95  55577       50277
## MSQ6_AG_517         0.99         0.99  67008      0.97  67008       63618
## MSQ6_AG_518         0.98         0.99  68679      0.95  68679       60000
## MSQ6_AG_519         0.99         0.99  59847      0.94  59847       49010
## MSQ6_AG_520         0.98         0.99  63359      0.93  63359       53096
## MSQ6_AG_521         0.95         0.99  19263      0.89  19263       15459
## MSQ6_AG_522         0.98         0.99  65289      0.95  65289       57063
## MSQ6_AG_523         0.98         0.99  58034      0.89  58034       43164
## MSQ6_AG_524         0.99         0.99  57373      0.91  57373       43837
## MSQ6_AG_525         0.99         0.99  72195      0.95  72195       56134
## MSQ6_AG_526         0.98         0.99  52673      0.88  52673       37386
## MSQ6_AG_527         0.99         0.99  60887      0.95  60887       48463
## MSQ6_AG_528         0.98         0.99  55620      0.91  55620       43056
## MSQ6_AG_529         0.97         0.99  35722      0.91  35722       27699
## MSQ6_AG_530         0.98         0.99  57328      0.86  57328       43663
##             length_filtered tabled_pc chimera_out_pc length_filtered_pc Sample
## MSQ4_AG_012           63132         1           0.75                  1  S_012
## MSQ4_AG_013           28004         1           0.78                  1  S_013
## MSQ4_AG_014           34919         1           0.73                  1  S_014
## MSQ4_AG_015           36310         1           0.72                  1  S_015
## MSQ4_AG_016           25798         1           0.75                  1  S_016
## MSQ4_AG_017           31556         1           0.71                  1  S_017
## MSQ4_AG_018           39875         1           0.88                  1  S_018
## MSQ4_AG_019           33557         1           0.76                  1  S_019
## MSQ4_AG_020           43716         1           0.70                  1  S_020
## MSQ4_AG_021           23176         1           0.79                  1  S_021
## MSQ4_AG_022           41438         1           0.83                  1  S_022
## MSQ4_AG_023           40845         1           0.84                  1  S_023
## MSQ4_AG_024           28168         1           0.70                  1  S_024
## MSQ4_AG_025           38257         1           0.74                  1  S_025
## MSQ4_AG_027           36624         1           0.69                  1  S_027
## MSQ4_AG_028           26811         1           0.72                  1  S_028
## MSQ4_AG_029           24267         1           0.70                  1  S_029
## MSQ4_AG_030           47010         1           0.75                  1  S_030
## MSQ4_AG_031           38050         1           0.78                  1  S_031
## MSQ4_AG_033           44462         1           0.69                  1  S_033
## MSQ4_AG_034           33723         1           0.72                  1  S_034
## MSQ4_AG_035           48208         1           0.81                  1  S_035
## MSQ4_AG_036           47982         1           0.71                  1  S_036
## MSQ4_AG_037           21250         1           0.71                  1  S_037
## MSQ4_AG_038           34480         1           0.74                  1  S_038
## MSQ4_AG_039           28166         1           0.69                  1  S_039
## MSQ4_AG_040           23637         1           0.74                  1  S_040
## MSQ4_AG_041           28883         1           0.73                  1  S_041
## MSQ4_AG_042           36332         1           0.77                  1  S_042
## MSQ4_AG_043           36253         1           0.77                  1  S_043
## MSQ4_AG_044           33676         1           0.74                  1  S_044
## MSQ4_AG_045           22118         1           0.67                  1  S_045
## MSQ4_AG_046           29882         1           0.71                  1  S_046
## MSQ4_AG_047           35244         1           0.77                  1  S_047
## MSQ4_AG_049           37365         1           0.71                  1  S_049
## MSQ4_AG_050           33113         1           0.74                  1  S_050
## MSQ4_AG_051           33451         1           0.69                  1  S_051
## MSQ4_AG_052           43344         1           0.76                  1  S_052
## MSQ4_AG_053           26692         1           0.70                  1  S_053
## MSQ4_AG_054           35882         1           0.67                  1  S_054
## MSQ4_AG_055           34369         1           0.68                  1  S_055
## MSQ4_AG_056           18872         1           0.68                  1  S_056
## MSQ4_AG_057           53231         1           0.89                  1  S_057
## MSQ4_AG_058           46387         1           0.83                  1  S_058
## MSQ4_AG_059           52456         1           0.91                  1  S_059
## MSQ4_AG_062           78943         1           0.93                  1  S_060
## MSQ4_AG_064           42321         1           0.94                  1  S_061
## MSQ4_AG_065           43474         1           0.92                  1  S_062
## MSQ4_AG_072           24746         1           0.95                  1  S_063
## MSQ4_AG_076           55585         1           0.95                  1  S_064
## MSQ4_AG_077           37756         1           0.88                  1  S_065
## MSQ4_AG_080           36761         1           0.95                  1  S_066
## MSQ4_AG_083           39045         1           0.77                  1  S_067
## MSQ4_AG_084           45763         1           0.77                  1  S_068
## MSQ4_AG_085           31589         1           0.74                  1  S_069
## MSQ4_AG_086           47584         1           0.77                  1  S_070
## MSQ4_AG_087           40651         1           0.88                  1  S_071
## MSQ4_AG_088           32150         1           0.87                  1  S_072
## MSQ4_AG_089           45583         1           0.84                  1  S_073
## MSQ4_AG_091           57555         1           1.00                  1  S_074
## MSQ4_AG_092           50090         1           0.79                  1  S_075
## MSQ4_AG_093           35237         1           0.76                  1  S_076
## MSQ4_AG_095           37534         1           0.78                  1  S_077
## MSQ4_AG_096           34866         1           0.65                  1  S_078
## MSQ4_AG_101           32686         1           0.83                  1  S_080
## MSQ4_AG_103           30485         1           0.78                  1  S_082
## MSQ4_AG_104           28858         1           0.76                  1  S_083
## MSQ4_AG_105           30798         1           0.75                  1  S_084
## MSQ4_AG_106           18531         1           0.76                  1  S_085
## MSQ4_AG_107           26957         1           0.73                  1  S_086
## MSQ4_AG_108           27091         1           0.77                  1  S_087
## MSQ4_AG_109           29029         1           0.77                  1  S_088
## MSQ4_AG_111           26044         1           0.71                  1  S_089
## MSQ4_AG_112           45944         1           0.92                  1  S_090
## MSQ4_AG_113           36643         1           0.85                  1  S_091
## MSQ4_AG_117           30639         1           0.74                  1  S_093
## MSQ4_AG_119           27533         1           0.78                  1  S_095
## MSQ4_AG_121           35845         1           0.72                  1  S_096
## MSQ4_AG_125           29264         1           0.81                  1  S_098
## MSQ4_AG_127           27454         1           0.77                  1  S_100
## MSQ4_AG_128           45962         1           0.80                  1  S_101
## MSQ4_AG_129           34436         1           0.73                  1  S_102
## MSQ4_AG_130           30300         1           0.76                  1  S_103
## MSQ4_AG_131           28664         1           0.68                  1  S_104
## MSQ4_AG_132           44270         1           0.78                  1  S_105
## MSQ4_AG_133           47432         1           0.92                  1  S_106
## MSQ4_AG_135           39379         1           0.88                  1  S_108
## MSQ4_AG_136           34867         1           0.76                  1  S_109
## MSQ4_AG_141           28476         1           0.77                  1  S_111
## MSQ4_AG_143           35175         1           0.76                  1  S_113
## MSQ4_AG_144           35801         1           0.69                  1  S_114
## MSQ4_AG_145           36833         1           0.67                  1  S_115
## MSQ4_AG_148           68895         1           0.73                  1  S_117
## MSQ4_AG_149           41929         1           0.74                  1  S_118
## MSQ4_AG_151           31350         1           0.73                  1  S_120
## MSQ4_AG_152           39031         1           0.73                  1  S_121
## MSQ4_AG_153           49753         1           0.75                  1  S_122
## MSQ4_AG_154           38218         1           0.74                  1  S_123
## MSQ4_AG_155           30522         1           0.67                  1  S_124
## MSQ4_AG_156           33495         1           0.74                  1  S_125
## MSQ4_AG_157           34703         1           0.76                  1  S_126
## MSQ4_AG_159           34623         1           0.78                  1  S_128
## MSQ4_AG_160           34466         1           0.74                  1  S_129
## MSQ4_AG_161           36003         1           0.76                  1  S_130
## MSQ4_AG_162           18381         1           0.75                  1  S_131
## MSQ4_AG_163           26906         1           0.75                  1  S_132
## MSQ4_AG_164           40076         1           0.76                  1  S_133
## MSQ4_AG_165           31663         1           0.78                  1  S_134
## MSQ4_AG_167           27235         1           0.76                  1  S_135
## MSQ4_AG_168           38492         1           0.72                  1  S_136
## MSQ4_AG_169           35599         1           0.76                  1  S_137
## MSQ4_AG_170           34791         1           0.74                  1  S_138
## MSQ4_AG_171           34171         1           0.70                  1  S_139
## MSQ4_AG_172           33548         1           0.75                  1  S_140
## MSQ4_AG_173           33771         1           0.72                  1  S_141
## MSQ4_AG_175           37348         1           0.66                  1  S_142
## MSQ4_AG_176           38771         1           0.80                  1  S_143
## MSQ4_AG_177           29662         1           0.74                  1  S_144
## MSQ4_AG_178           24028         1           0.68                  1  S_145
## MSQ4_AG_179           27466         1           0.75                  1  S_146
## MSQ4_AG_180           31624         1           0.80                  1  S_147
## MSQ4_AG_181           28373         1           0.78                  1  S_148
## MSQ4_AG_183           28511         1           0.77                  1  S_150
## MSQ4_AG_184           30442         1           0.71                  1  S_151
## MSQ4_AG_185           31339         1           0.74                  1  S_152
## MSQ4_AG_186           31727         1           0.75                  1  S_153
## MSQ4_AG_187           23085         1           0.74                  1  S_154
## MSQ4_AG_188           25966         1           0.76                  1  S_155
## MSQ4_AG_189           34953         1           0.72                  1  S_156
## MSQ4_AG_191           38217         1           0.74                  1  S_158
## MSQ4_AG_192           37249         1           0.85                  1  S_159
## MSQ4_AG_193           39378         1           0.88                  1  S_160
## MSQ4_AG_194           32588         1           0.79                  1  S_161
## MSQ4_AG_195           24980         1           0.75                  1  S_162
## MSQ4_AG_196           32211         1           0.73                  1  S_163
## MSQ4_AG_197           31808         1           0.74                  1  S_164
## MSQ4_AG_198           21792         1           0.79                  1  S_165
## MSQ4_AG_201           12911         1           0.76                  1  S_168
## MSQ4_AG_205           22030         1           0.88                  1  S_172
## MSQ4_AG_207           35316         1           0.80                  1  S_174
## MSQ4_AG_208           31211         1           0.71                  1  S_175
## MSQ4_AG_209           31033         1           0.72                  1  S_176
## MSQ4_AG_210           32410         1           0.79                  1  S_177
## MSQ4_AG_211           22976         1           0.76                  1  S_178
## MSQ4_AG_212           28577         1           0.76                  1  S_179
## MSQ4_AG_213           30989         1           0.77                  1  S_180
## MSQ4_AG_215           29770         1           0.66                  1  S_182
## MSQ4_AG_216           31999         1           0.70                  1  S_183
## MSQ4_AG_217           25845         1           0.62                  1  S_184
## MSQ4_AG_218           29813         1           0.63                  1  S_185
## MSQ4_AG_219           29063         1           0.78                  1  S_186
## MSQ4_AG_220           27283         1           0.66                  1  S_187
## MSQ4_AG_221           28554         1           0.67                  1  S_188
## MSQ4_AG_223          121734         1           0.73                  1  S_190
## MSQ4_AG_224           43308         1           0.76                  1  S_191
## MSQ4_AG_225           44919         1           0.74                  1  S_192
## MSQ4_AG_226           42654         1           0.73                  1  S_193
## MSQ4_AG_227           30151         1           0.68                  1  S_194
## MSQ4_AG_228           26182         1           0.74                  1  S_195
## MSQ4_AG_229           34084         1           0.78                  1  S_196
## MSQ4_AG_230           27509         1           0.66                  1  S_197
## MSQ4_AG_231           29257         1           0.73                  1  S_198
## MSQ4_AG_232           39903         1           0.89                  1  S_199
## MSQ4_AG_233           31217         1           0.72                  1  S_200
## MSQ4_AG_234           36152         1           0.71                  1  S_201
## MSQ4_AG_235           34859         1           0.78                  1  S_202
## MSQ4_AG_236           35312         1           0.82                  1  S_203
## MSQ4_AG_237           27057         1           0.71                  1  S_204
## MSQ4_AG_238           25671         1           0.79                  1  S_205
## MSQ4_AG_239           38863         1           0.81                  1  S_206
## MSQ4_AG_240           37376         1           0.75                  1  S_207
## MSQ4_AG_241           28959         1           0.65                  1  S_208
## MSQ4_AG_242           36918         1           0.68                  1  S_209
## MSQ4_AG_243           22034         1           0.67                  1  S_210
## MSQ4_AG_244           35309         1           0.72                  1  S_211
## MSQ4_AG_245           25726         1           0.63                  1  S_212
## MSQ4_AG_247           39348         1           0.80                  1  S_214
## MSQ4_AG_248           33598         1           0.67                  1  S_215
## MSQ4_AG_249           33624         1           0.75                  1  S_216
## MSQ6_AG_250           26929         1           0.71                  1  S_219
## MSQ6_AG_251           44125         1           0.69                  1  S_220
## MSQ6_AG_252           49814         1           0.73                  1  S_221
## MSQ6_AG_253           42532         1           0.77                  1  S_222
## MSQ6_AG_254           49859         1           0.73                  1  S_223
## MSQ6_AG_255           44187         1           0.73                  1  S_224
## MSQ6_AG_256          124781         1           0.72                  1  S_225
## MSQ6_AG_257           39110         1           0.73                  1  S_226
## MSQ6_AG_258           47209         1           0.72                  1  S_227
## MSQ6_AG_259           39238         1           0.75                  1  S_228
## MSQ6_AG_260           49569         1           0.83                  1  S_229
## MSQ6_AG_261           32999         1           0.73                  1  S_230
## MSQ6_AG_262           89836         1           0.79                  1  S_231
## MSQ6_AG_263           37806         1           0.70                  1  S_232
## MSQ6_AG_264           44726         1           0.76                  1  S_233
## MSQ6_AG_265           37761         1           0.72                  1  S_234
## MSQ6_AG_266           32968         1           0.72                  1  S_235
## MSQ6_AG_267           29158         1           0.71                  1  S_236
## MSQ6_AG_268           37290         1           0.73                  1  S_237
## MSQ6_AG_269           24804         1           0.73                  1  S_238
## MSQ6_AG_270           45108         1           0.73                  1  S_239
## MSQ6_AG_271           35201         1           0.67                  1  S_240
## MSQ6_AG_272           37456         1           0.73                  1  S_241
## MSQ6_AG_273           33862         1           0.68                  1  S_242
## MSQ6_AG_274           18014         1           0.74                  1  S_243
## MSQ6_AG_275           43190         1           0.72                  1  S_244
## MSQ6_AG_276           45910         1           0.70                  1  S_245
## MSQ6_AG_277           36000         1           0.72                  1  S_246
## MSQ6_AG_278           34478         1           0.73                  1  S_247
## MSQ6_AG_279           36188         1           0.70                  1  S_248
## MSQ6_AG_281           50265         1           0.74                  1  S_250
## MSQ6_AG_282           36121         1           0.85                  1  S_251
## MSQ6_AG_283           27664         1           0.74                  1  S_252
## MSQ6_AG_284           39123         1           0.70                  1  S_253
## MSQ6_AG_285           34414         1           0.72                  1  S_254
## MSQ6_AG_286           41213         1           0.74                  1  S_255
## MSQ6_AG_287           39884         1           0.75                  1  S_256
## MSQ6_AG_288           35603         1           0.74                  1  S_257
## MSQ6_AG_290           22393         1           0.76                  1  S_259
## MSQ6_AG_291           37401         1           0.70                  1  S_260
## MSQ6_AG_292           44075         1           0.69                  1  S_261
## MSQ6_AG_293           33102         1           0.77                  1  S_262
## MSQ6_AG_294           39177         1           0.74                  1  S_263
## MSQ6_AG_295           36371         1           0.71                  1  S_264
## MSQ6_AG_296           33074         1           0.70                  1  S_265
## MSQ6_AG_298           25918         1           0.78                  1  S_267
## MSQ6_AG_299           44286         1           0.75                  1  S_268
## MSQ6_AG_300           41040         1           0.74                  1  S_269
## MSQ6_AG_301           33787         1           0.72                  1  S_270
## MSQ6_AG_302           49311         1           0.79                  1  S_271
## MSQ6_AG_303           47625         1           0.75                  1  S_272
## MSQ6_AG_304           38035         1           0.78                  1  S_273
## MSQ6_AG_305           59937         1           0.74                  1  S_274
## MSQ6_AG_306           30890         1           0.73                  1  S_275
## MSQ6_AG_307           36986         1           0.78                  1  S_276
## MSQ6_AG_308           40691         1           0.77                  1  S_277
## MSQ6_AG_309           31897         1           0.77                  1  S_278
## MSQ6_AG_310           47728         1           0.75                  1  S_279
## MSQ6_AG_311           37306         1           0.79                  1  S_280
## MSQ6_AG_312           32803         1           0.76                  1  S_281
## MSQ6_AG_313           39652         1           0.74                  1  S_282
## MSQ6_AG_314           39270         1           0.75                  1  S_283
## MSQ6_AG_315           46336         1           0.77                  1  S_284
## MSQ6_AG_316           45015         1           0.86                  1  S_285
## MSQ6_AG_317           30947         1           0.76                  1  S_286
## MSQ6_AG_318           46795         1           0.75                  1  S_287
## MSQ6_AG_319           43420         1           0.72                  1  S_288
## MSQ6_AG_320           36090         1           0.79                  1  S_289
## MSQ6_AG_321           49278         1           0.80                  1  S_290
## MSQ6_AG_322           30670         1           0.77                  1  S_291
## MSQ6_AG_323           36889         1           0.78                  1  S_292
## MSQ6_AG_324           40654         1           0.84                  1  S_293
## MSQ6_AG_325           35242         1           0.79                  1  S_294
## MSQ6_AG_326           44122         1           0.83                  1  S_295
## MSQ6_AG_327           35011         1           0.78                  1  S_296
## MSQ6_AG_328           36199         1           0.76                  1  S_297
## MSQ6_AG_329           32872         1           0.80                  1  S_298
## MSQ6_AG_330           25123         1           0.76                  1  S_299
## MSQ6_AG_331           34215         1           0.77                  1  S_300
## MSQ6_AG_332           39342         1           0.77                  1  S_301
## MSQ6_AG_333           27584         1           0.76                  1  S_302
## MSQ6_AG_334           36753         1           0.71                  1  S_303
## MSQ6_AG_335           34090         1           0.67                  1  S_304
## MSQ6_AG_336           32369         1           0.74                  1  S_305
## MSQ6_AG_337           33486         1           0.72                  1  S_306
## MSQ6_AG_339           37026         1           0.73                  1  S_308
## MSQ6_AG_340           30284         1           0.70                  1  S_309
## MSQ6_AG_341           24770         1           0.74                  1  S_310
## MSQ6_AG_342           41706         1           0.71                  1  S_311
## MSQ6_AG_343           34794         1           0.74                  1  S_312
## MSQ6_AG_344           38160         1           0.79                  1  S_313
## MSQ6_AG_346           37273         1           0.74                  1  S_315
## MSQ6_AG_347           35747         1           0.74                  1  S_316
## MSQ6_AG_348           36364         1           0.72                  1  S_317
## MSQ6_AG_349           32689         1           0.76                  1  S_318
## MSQ6_AG_350           35016         1           0.75                  1  S_319
## MSQ6_AG_351           33891         1           0.74                  1  S_320
## MSQ6_AG_352           35926         1           0.73                  1  S_321
## MSQ6_AG_353           33868         1           0.72                  1  S_322
## MSQ6_AG_354           38224         1           0.76                  1  S_323
## MSQ6_AG_355           41368         1           0.83                  1  S_324
## MSQ6_AG_356           29213         1           0.74                  1  S_325
## MSQ6_AG_357           37943         1           0.76                  1  S_326
## MSQ6_AG_358           33353         1           0.81                  1  S_327
## MSQ6_AG_359           41571         1           0.75                  1  S_328
## MSQ6_AG_360           39599         1           0.74                  1  S_329
## MSQ6_AG_361           32148         1           0.77                  1  S_330
## MSQ6_AG_362           35602         1           0.77                  1  S_331
## MSQ6_AG_363           27164         1           0.77                  1  S_332
## MSQ6_AG_364           35941         1           0.82                  1  S_333
## MSQ6_AG_365           30861         1           0.79                  1  S_334
## MSQ6_AG_366           37278         1           0.80                  1  S_335
## MSQ6_AG_367           34694         1           0.80                  1  S_336
## MSQ6_AG_368           37761         1           0.76                  1  S_337
## MSQ6_AG_369           38432         1           0.83                  1  S_338
## MSQ6_AG_370           35563         1           0.77                  1  S_339
## MSQ6_AG_371           35425         1           0.79                  1  S_340
## MSQ6_AG_372           37388         1           0.78                  1  S_341
## MSQ6_AG_373           30253         1           0.84                  1  S_342
## MSQ6_AG_374           22578         1           0.79                  1  S_343
## MSQ6_AG_375           41121         1           0.79                  1  S_344
## MSQ6_AG_376           18120         1           0.76                  1  S_345
## MSQ6_AG_377           30174         1           0.78                  1  S_346
## MSQ6_AG_378           34059         1           0.74                  1  S_347
## MSQ6_AG_379           36536         1           0.78                  1  S_348
## MSQ6_AG_380           54422         1           0.72                  1  S_349
## MSQ6_AG_381           33745         1           0.70                  1  S_350
## MSQ6_AG_382           47662         1           0.70                  1  S_351
## MSQ6_AG_383           45958         1           0.73                  1  S_352
## MSQ6_AG_384           45184         1           0.80                  1  S_353
## MSQ6_AG_385           34946         1           0.73                  1  S_354
## MSQ6_AG_386           32861         1           0.71                  1  S_355
## MSQ6_AG_387           29841         1           0.72                  1  S_356
## MSQ6_AG_388           32762         1           0.76                  1  S_357
## MSQ6_AG_389           42292         1           0.78                  1  S_358
## MSQ6_AG_390           34330         1           0.76                  1  S_359
## MSQ6_AG_391           45537         1           0.78                  1  S_360
## MSQ6_AG_392           34869         1           0.77                  1  S_361
## MSQ6_AG_393           34710         1           0.79                  1  S_362
## MSQ6_AG_394           42780         1           0.77                  1  S_363
## MSQ6_AG_395           35991         1           0.75                  1  S_364
## MSQ6_AG_396           45382         1           0.78                  1  S_365
## MSQ6_AG_397           34890         1           0.80                  1  S_366
## MSQ6_AG_398           32099         1           0.76                  1  S_367
## MSQ6_AG_399           60535         1           0.76                  1  S_368
## MSQ6_AG_400           49432         1           0.76                  1  S_369
## MSQ6_AG_401           41237         1           0.75                  1  S_370
## MSQ6_AG_402           35908         1           0.78                  1  S_371
## MSQ6_AG_403           29153         1           0.74                  1  S_372
## MSQ6_AG_404           31418         1           0.75                  1  S_373
## MSQ6_AG_405           35131         1           0.77                  1  S_374
## MSQ6_AG_406           41200         1           0.79                  1  S_375
## MSQ6_AG_407           61652         1           0.78                  1  S_376
## MSQ6_AG_408           51398         1           0.77                  1  S_377
## MSQ6_AG_409           36318         1           0.78                  1  S_378
## MSQ6_AG_410           31866         1           0.79                  1  S_379
## MSQ6_AG_411           24987         1           0.75                  1  S_380
## MSQ6_AG_412           40913         1           0.78                  1  S_381
## MSQ6_AG_413           33713         1           0.82                  1  S_382
## MSQ6_AG_414           31111         1           0.77                  1  S_383
## MSQ6_AG_415           42346         1           0.77                  1  S_384
## MSQ6_AG_416           37759         1           0.77                  1  S_385
## MSQ6_AG_417           28740         1           0.77                  1  S_386
## MSQ6_AG_418           36412         1           0.77                  1  S_387
## MSQ6_AG_419           30021         1           0.73                  1  S_388
## MSQ6_AG_420           32891         1           0.78                  1  S_389
## MSQ6_AG_421           32405         1           0.77                  1  S_390
## MSQ6_AG_422           40869         1           0.80                  1  S_391
## MSQ6_AG_423           33271         1           0.74                  1  S_392
## MSQ6_AG_424           29615         1           0.76                  1  S_393
## MSQ6_AG_425           31675         1           0.75                  1  S_394
## MSQ6_AG_426           37581         1           0.74                  1  S_395
## MSQ6_AG_427           36710         1           0.76                  1  S_396
## MSQ6_AG_428           39430         1           0.75                  1  S_397
## MSQ6_AG_429           40730         1           0.79                  1  S_398
## MSQ6_AG_430           34324         1           0.76                  1  S_399
## MSQ6_AG_431           47958         1           0.86                  1  S_400
## MSQ6_AG_432           42098         1           0.79                  1  S_401
## MSQ6_AG_433           36215         1           0.77                  1  S_402
## MSQ6_AG_434           30546         1           0.75                  1  S_403
## MSQ6_AG_490           52584         1           0.87                  1  S_438
## MSQ6_AG_491           50751         1           0.84                  1  S_439
## MSQ6_AG_492           33283         1           0.86                  1  S_440
## MSQ6_AG_493           58699         1           0.88                  1  S_441
## MSQ6_AG_494           58881         1           0.98                  1  S_442
## MSQ6_AG_495           50758         1           0.83                  1  S_443
## MSQ6_AG_496           60863         1           0.96                  1  S_444
## MSQ6_AG_497           36375         1           0.85                  1  S_445
## MSQ6_AG_498           50918         1           0.83                  1  S_446
## MSQ6_AG_499           58042         1           0.92                  1  S_447
## MSQ6_AG_500           46975         1           0.87                  1  S_448
## MSQ6_AG_501           57531         1           0.87                  1  S_449
## MSQ6_AG_502           53894         1           0.87                  1  S_450
## MSQ6_AG_503           43444         1           0.99                  1  S_451
## MSQ6_AG_504           63433         1           0.86                  1  S_452
## MSQ6_AG_505           14674         1           0.86                  1  S_453
## MSQ6_AG_506           19782         1           0.89                  1  S_454
## MSQ6_AG_509           15767         1           0.86                  1  S_457
## MSQ6_AG_511           22168         1           0.87                  1  S_459
## MSQ6_AG_512           22917         1           0.86                  1  S_460
## MSQ6_AG_514           56620         1           0.82                  1  S_462
## MSQ6_AG_515           68258         1           0.91                  1  S_463
## MSQ6_AG_516           50277         1           0.90                  1  S_464
## MSQ6_AG_517           63618         1           0.95                  1  S_465
## MSQ6_AG_518           60000         1           0.87                  1  S_466
## MSQ6_AG_519           49010         1           0.82                  1  S_467
## MSQ6_AG_520           53096         1           0.84                  1  S_468
## MSQ6_AG_521           15459         1           0.80                  1  S_469
## MSQ6_AG_522           57063         1           0.87                  1  S_470
## MSQ6_AG_523           43164         1           0.74                  1  S_471
## MSQ6_AG_524           43837         1           0.76                  1  S_472
## MSQ6_AG_525           56134         1           0.78                  1  S_473
## MSQ6_AG_526           37384         1           0.71                  1  S_474
## MSQ6_AG_527           48463         1           0.80                  1  S_475
## MSQ6_AG_528           43056         1           0.77                  1  S_476
## MSQ6_AG_529           27699         1           0.78                  1  S_477
## MSQ6_AG_530           43663         1           0.76                  1  S_478
##              reads remove responsible_seq owner region_1  run experiment
## MSQ4_AG_012  48960      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_013  20385      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_014  23342      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_015  23590      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_016  18212      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_017  22353      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_018  31561      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_019  25224      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_020  28421      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_021  16768      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_022  27781      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_023  26325      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_024  22660      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_025  25542      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_027  31281      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_028  13506      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_029  13524      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_030  41478      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_031  27761      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_033  29209      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_034  22710      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_035  37399      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_036  29296      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_037  18759      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_038  29865      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_039  28061      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_040  16779      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_041  18461      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_042  28067      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_043  24041      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_044  32155      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_045  18910      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_046  26743      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_047  21306      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_049  22418      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_050  27441      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_051  30582      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_052  31698      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_053  17888      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_054  26929      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_055  29620      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_056  10306      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_057  38560      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_058  31142      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_059  38293      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_062  46431      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_064  27225      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_065  26286      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_072  16896      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_076  34565      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_077  20024      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_080  20935      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_083  13463      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_084  16582      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_085  21248      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_086  29225      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_087  30111      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_088  25506      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_089  37671      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_091  42171      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_092  42724      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_093  26838      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_095  29438      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_096  30533      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_101  19850      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_103  21225      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_104  16463      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_105  22684      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_106  14078      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_107  18019      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_108  12907      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_109  18343      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_111  16565      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_112  29007      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_113  18743      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_117  17886      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_119  16730      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_121  16882      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_125  17420      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_127  22454      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_128  26337      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_129  30142      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_130  18481      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_131  18582      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_132  34709      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_133  41146      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_135  28965      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_136  27319      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_141  15448      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_143  20556      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_144  18792      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_145  21277      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_148  38364      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_149  32066      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_151  29288      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_152  22915      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_153  37904      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_154  30432      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_155  24624      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_156  26484      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_157  23585      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_159  26585      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_160  29184      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_161  27148      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_162  12458      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_163  24549      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_164  29877      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_165  17916      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_167  21697      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_168  36177      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_169  27337      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_170  26173      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_171  34590      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_172  29459      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_173  32652      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_175  37357      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_176  14564      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_177  23021      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_178  17451      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_179  17820      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_180  27594      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_181  18430      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_183  14229      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_184  20532      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_185  23171      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_186  22289      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_187  18409      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_188  21405      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_189  27654      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_191  27438      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_192  28350      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_193  31426      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_194  27764      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_195  18918      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_196  29009      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_197  27959      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_198  15314      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_201  11671      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_205  17452      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_207  28609      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_208  23553      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_209  24203      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_210  26296      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_211  17438      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_212  22574      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_213  26739      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_215  18875      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_216  23681      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_217  18749      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_218  18912      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_219  12735      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_220  18751      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_221  18213      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_223 105501      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_224  37502      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_225  40844      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_226  38679      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_227  21736      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_228  19344      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_229  25190      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_230  22205      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_231  22528      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_232  26132      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_233  26777      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_234  32832      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_235  26317      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_236  29962      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_237  21919      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_238  19521      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_239  28966      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_240  31434      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_241  17971      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_242  27908      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_243  11430      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_244  22644      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_245  13051      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_247  27635      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_248  19923      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ4_AG_249  26099      s        Annelies    AG       V4 MSQ4     mIMT_1
## MSQ6_AG_250  18616      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_251  27311      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_252  28127      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_253  23855      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_254  40793      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_255  37037      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_256  84206      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_257  24313      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_258  30945      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_259  17347      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_260  41233      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_261  17757      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_262  76765      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_263  24326      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_264  30767      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_265  23093      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_266  18061      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_267  17456      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_268  25508      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_269  11644      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_270  23015      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_271  21122      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_272  25009      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_273  18619      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_274  13480      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_275  25204      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_276  31518      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_277  23000      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_278  15369      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_279  24112      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_281  26861      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_282  31371      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_283  19402      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_284  26676      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_285  23704      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_286  29094      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_287  27534      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_288  22057      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_290  15048      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_291  25086      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_292  30359      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_293  22352      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_294  34683      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_295  22771      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_296  26224      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_298  19363      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_299  33058      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_300  30595      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_301  25220      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_302  39359      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_303  34494      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_304  29387      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_305  48529      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_306  20899      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_307  26411      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_308  33269      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_309  24903      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_310  42828      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_311  24629      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_312  28530      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_313  32978      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_314  33285      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_315  38605      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_316  35003      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_317  25237      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_318  43348      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_319  40817      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_320  30025      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_321  28521      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_322  19187      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_323  25485      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_324  32392      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_325  16445      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_326  31454      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_327  18757      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_328  19185      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_329  23689      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_330  18327      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_331  23353      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_332  27068      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_333  16367      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_334  34241      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_335  32389      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_336  23822      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_337  20350      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_339  31462      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_340  14098      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_341  20222      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_342  32489      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_343  21782      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_344  24703      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_346  29733      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_347  25822      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_348  25042      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_349  22374      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_350  21462      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_351  24471      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_352  25050      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_353  16390      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_354  26938      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_355  24248      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_356  17733      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_357  34400      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_358  23792      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_359  38358      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_360  32322      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_361  27549      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_362  26252      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_363  19361      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_364  24551      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_365  27389      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_366  29567      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_367  26725      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_368  31468      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_369  32157      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_370  23974      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_371  22451      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_372  26245      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_373  23351      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_374  16927      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_375  36013      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_376  16695      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_377  25765      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_378  24392      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_379  35828      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_380  42252      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_381  22422      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_382  35537      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_383  40142      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_384  39030      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_385  26885      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_386  23918      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_387  17693      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_388  23570      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_389  28400      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_390  31764      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_391  36005      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_392  30819      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_393  22694      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_394  29687      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_395  30662      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_396  35373      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_397  26152      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_398  29306      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_399  44910      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_400  33218      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_401  36986      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_402  26873      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_403  25001      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_404  26432      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_405  27462      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_406  30513      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_407  44225      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_408  35772      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_409  26830      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_410  24015      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_411  23950      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_412  34170      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_413  26234      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_414  25249      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_415  36547      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_416  32805      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_417  23332      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_418  27881      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_419  24558      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_420  29044      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_421  26611      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_422  35091      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_423  31460      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_424  24200      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_425  25743      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_426  37239      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_427  31120      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_428  31620      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_429  33578      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_430  28910      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_431  42838      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_432  36424      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_433  30466      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_434  21689      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_490  42522      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_491  39448      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_492  25183      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_493  41289      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_494  45501      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_495  32594      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_496  47040      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_497  28263      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_498  39328      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_499  48814      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_500  34845      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_501  53167      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_502  46455      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_503  36462      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_504  50834      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_505  11269      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_506  15473      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_509  11622      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_511  15183      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_512  17188      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_514  47195      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_515  56791      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_516  37290      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_517  54378      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_518  46183      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_519  40486      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_520  44884      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_521  12799      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_522  47385      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_523  38843      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_524  30577      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_525  33737      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_526  36823      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_527  31171      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_528  32578      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_529  17380      s        Annelies    AG       V4 MSQ6     mIMT_2
## MSQ6_AG_530  37773      s        Annelies    AG       V4 MSQ6     mIMT_2
##             Fermentation facility    day sample_type reactor day_IR
## MSQ4_AG_012         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_013         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_014         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_015         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_016         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_017         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_018         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_019         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_020         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_021         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_022         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_023         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_024         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_025         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_027         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_028         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_029         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_030         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_031         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_033         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_034         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_035         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_036         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_037         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_038         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_039         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_040         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_041         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_042         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_043         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_044         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_045         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_046         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ4_AG_047         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_049         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_050         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_051         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_052         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_053         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_054         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_055         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_056         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_057         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_058         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_059         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_062         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_064         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_065         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_072         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_076         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_077         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_080         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ4_AG_083         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_084         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_085         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_086         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_087         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_088         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_089         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_091         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_092         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_093         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_095         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_096         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_101         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_103         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_104         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_105         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_106         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_107         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_108         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_109         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_111         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_112         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_113         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ4_AG_117         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_119         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_121         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_125         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_127         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_128         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_129         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_130         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_131         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_132         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_133         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_135         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_136         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_141         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_143         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_144         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_145         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_148         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_149         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ4_AG_151         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_152         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_153         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_154         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_155         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_156         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_157         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_159         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_160         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_161         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_162         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_163         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_164         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_165         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_167         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_168         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_169         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_170         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_171         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_172         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_173         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_175         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_176         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_177         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_178         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_179         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_180         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_181         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_183         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ4_AG_184         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_185         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_186         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_187         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_188         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_189         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_191         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_192         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_193         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_194         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_195         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_196         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_197         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_198         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_201         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_205         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_207         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_208         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_209         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_210         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_211         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_212         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_213         <NA>   Borsig Day_14       Cecum    <NA>     NA
## MSQ4_AG_215         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_216         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_217         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_218         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_219         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_220         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_221         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_223         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_224         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_225         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_226         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_227         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_228         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_229         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_230         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_231         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_232         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_233         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_234         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_235         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_236         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_237         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_238         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_239         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_240         <NA>   Borsig Day_12       Feces    <NA>     NA
## MSQ4_AG_241         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ4_AG_242         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ4_AG_243         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ4_AG_244         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ4_AG_245         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ4_AG_247         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ4_AG_248         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ4_AG_249         <NA>   Borsig Day_-1       Feces    <NA>     NA
## MSQ6_AG_250         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_251         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_252         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_253         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_254         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_255         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_256         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_257         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_258         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_259         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_260         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_261         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_262         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_263         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_264         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_265         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_266         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_267         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_268         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_269         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_270         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_271         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_272         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_273         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_274         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_275         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_276         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_277         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_278         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_279         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_281         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_282         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_283         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_284         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_285         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_286         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_287         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_288         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_290         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_291         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_292         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_293         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_294         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_295         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_296         <NA>   Scharl Day_-1       Feces    <NA>     NA
## MSQ6_AG_298         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_299         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_300         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_301         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_302         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_303         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_304         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_305         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_306         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_307         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_308         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_309         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_310         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_311         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_312         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_313         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_314         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_315         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_316         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_317         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_318         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_319         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_320         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_321         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_322         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_323         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_324         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_325         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_326         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_327         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_328         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_329         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_330         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_331         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_332         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_333         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_334         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_335         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_336         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_337         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_339         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_340         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_341         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_342         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_343         <NA>   Scharl Day_12       Feces    <NA>     NA
## MSQ6_AG_344         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_346         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_347         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_348         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_349         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_350         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_351         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_352         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_353         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_354         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_355         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_356         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_357         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_358         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_359         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_360         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_361         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_362         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_363         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_364         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_365         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_366         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_367         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_368         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_369         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_370         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_371         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_372         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_373         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_374         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_375         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_376         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_377         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_378         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_379         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_380         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_381         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_382         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_383         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_384         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_385         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_386         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_387         <NA>   Scharl Day_14       Feces    <NA>     NA
## MSQ6_AG_388         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_389         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_390         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_391         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_392         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_393         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_394         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_395         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_396         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_397         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_398         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_399         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_400         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_401         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_402         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_403         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_404         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_405         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_406         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_407         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_408         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_409         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_410         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_411         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_412         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_413         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_414         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_415         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_416         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_417         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_418         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_419         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_420         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_421         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_422         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_423         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_424         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_425         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_426         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_427         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_428         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_429         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_430         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_431         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_432         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_433         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_434         <NA>   Scharl Day_14       Cecum    <NA>     NA
## MSQ6_AG_490         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_491         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_492         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_493         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_494         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_495         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_496         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_497         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_498         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_499         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_500         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_501         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_502         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_503         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_504         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_505         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_506         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_509         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_511         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_512         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_514         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_515         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_516         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_517         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_518         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_519         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_520         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_521         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_522         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_523         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_524         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_525         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_526         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_527         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_528         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_529         <NA>   Scharl  Day_6       Feces    <NA>     NA
## MSQ6_AG_530         <NA>   Scharl  Day_6       Feces    <NA>     NA
##             day_fermentation treatment_period sample_label treatment_exvivo
## MSQ4_AG_012             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_013             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_014             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_015             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_016             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_017             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_018             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_019             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_020             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_021             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_022             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_023             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_024             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_025             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_027             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_028             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_029             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_030             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_031             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_033             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_034             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_035             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_036             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_037             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_038             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_039             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_040             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_041             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_042             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_043             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_044             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_045             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_046             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_047             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_049             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_050             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_051             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_052             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_053             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_054             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_055             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_056             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_057             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_058             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_059             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_062             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_064             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_065             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_072             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_076             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_077             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_080             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_083             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_084             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_085             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_086             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_087             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_088             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_089             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_091             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_092             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_093             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_095             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_096             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_101             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_103             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_104             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_105             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_106             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_107             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_108             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_109             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_111             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_112             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_113             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_117             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_119             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_121             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_125             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_127             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_128             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_129             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_130             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_131             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_132             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_133             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_135             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_136             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_141             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_143             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_144             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_145             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_148             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_149             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_151             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_152             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_153             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_154             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_155             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_156             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_157             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_159             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_160             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_161             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_162             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_163             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_164             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_165             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_167             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_168             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_169             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_170             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_171             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_172             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_173             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_175             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_176             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_177             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_178             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_179             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_180             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_181             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_183             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_184             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_185             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_186             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_187             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_188             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_189             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_191             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_192             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_193             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_194             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_195             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_196             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_197             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_198             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_201             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_205             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_207             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_208             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_209             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_210             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_211             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_212             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_213             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_215             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_216             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_217             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_218             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_219             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_220             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_221             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_223             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_224             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_225             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_226             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_227             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_228             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_229             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_230             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_231             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_232             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_233             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_234             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_235             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_236             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_237             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_238             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_239             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_240             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_241             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_242             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_243             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_244             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_245             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_247             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_248             <NA>             <NA>         <NA>             <NA>
## MSQ4_AG_249             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_250             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_251             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_252             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_253             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_254             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_255             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_256             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_257             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_258             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_259             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_260             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_261             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_262             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_263             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_264             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_265             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_266             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_267             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_268             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_269             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_270             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_271             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_272             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_273             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_274             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_275             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_276             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_277             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_278             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_279             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_281             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_282             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_283             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_284             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_285             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_286             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_287             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_288             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_290             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_291             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_292             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_293             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_294             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_295             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_296             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_298             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_299             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_300             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_301             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_302             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_303             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_304             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_305             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_306             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_307             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_308             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_309             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_310             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_311             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_312             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_313             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_314             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_315             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_316             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_317             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_318             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_319             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_320             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_321             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_322             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_323             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_324             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_325             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_326             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_327             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_328             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_329             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_330             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_331             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_332             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_333             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_334             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_335             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_336             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_337             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_339             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_340             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_341             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_342             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_343             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_344             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_346             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_347             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_348             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_349             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_350             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_351             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_352             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_353             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_354             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_355             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_356             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_357             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_358             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_359             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_360             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_361             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_362             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_363             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_364             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_365             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_366             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_367             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_368             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_369             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_370             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_371             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_372             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_373             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_374             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_375             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_376             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_377             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_378             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_379             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_380             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_381             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_382             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_383             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_384             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_385             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_386             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_387             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_388             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_389             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_390             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_391             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_392             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_393             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_394             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_395             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_396             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_397             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_398             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_399             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_400             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_401             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_402             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_403             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_404             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_405             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_406             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_407             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_408             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_409             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_410             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_411             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_412             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_413             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_414             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_415             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_416             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_417             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_418             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_419             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_420             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_421             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_422             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_423             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_424             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_425             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_426             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_427             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_428             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_429             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_430             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_431             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_432             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_433             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_434             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_490             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_491             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_492             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_493             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_494             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_495             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_496             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_497             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_498             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_499             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_500             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_501             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_502             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_503             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_504             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_505             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_506             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_509             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_511             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_512             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_514             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_515             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_516             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_517             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_518             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_519             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_520             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_521             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_522             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_523             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_524             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_525             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_526             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_527             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_528             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_529             <NA>             <NA>         <NA>             <NA>
## MSQ6_AG_530             <NA>             <NA>         <NA>             <NA>
##               treatment mouse_label    BW    BW_percent BW_delta DAI
## MSQ4_AG_012         H2O         m_6 20.84   0.000000000     0.00  NA
## MSQ4_AG_013         H2O        m_19 23.99   0.000000000     0.00  NA
## MSQ4_AG_014         H2O        m_31 25.35   0.000000000     0.00  NA
## MSQ4_AG_015         H2O        m_33 27.54   0.000000000     0.00  NA
## MSQ4_AG_016         H2O        m_34 27.19   0.000000000     0.00  NA
## MSQ4_AG_017      H2O_co        m_10 21.24   0.000000000     0.00  NA
## MSQ4_AG_018      H2O_co        m_26 23.14   0.000000000     0.00  NA
## MSQ4_AG_019      H2O_co        m_27 23.09   0.000000000     0.00  NA
## MSQ4_AG_020      H2O_co        m_29 22.59   0.000000000     0.00  NA
## MSQ4_AG_021      H2O_co        m_30 22.00   0.000000000     0.00  NA
## MSQ4_AG_022         DSS         m_4 20.38   0.000000000     0.00  NA
## MSQ4_AG_023         DSS        m_16 22.00   0.000000000     0.00  NA
## MSQ4_AG_024         DSS        m_17 22.99   0.000000000     0.00  NA
## MSQ4_AG_025         DSS        m_18 23.56   0.000000000     0.00  NA
## MSQ4_AG_027      DSS_co        m_11 23.96   0.000000000     0.00  NA
## MSQ4_AG_028      DSS_co        m_21 20.22   0.000000000     0.00  NA
## MSQ4_AG_029      DSS_co        m_22 22.47   0.000000000     0.00  NA
## MSQ4_AG_030      DSS_co        m_23 23.47   0.000000000     0.00  NA
## MSQ4_AG_031      DSS_co        m_32 22.97   0.000000000     0.00  NA
## MSQ4_AG_033    Eubiotic         m_2 22.45   0.000000000     0.00  NA
## MSQ4_AG_034    Eubiotic         m_3 22.25   0.000000000     0.00  NA
## MSQ4_AG_035    Eubiotic         m_5 22.43   0.000000000     0.00  NA
## MSQ4_AG_036    Eubiotic        m_35 22.59   0.000000000     0.00  NA
## MSQ4_AG_037   Dysbiotic         m_7 21.13   0.000000000     0.00  NA
## MSQ4_AG_038   Dysbiotic         m_8 23.39   0.000000000     0.00  NA
## MSQ4_AG_039   Dysbiotic         m_9 22.84   0.000000000     0.00  NA
## MSQ4_AG_040   Dysbiotic        m_24 22.45   0.000000000     0.00  NA
## MSQ4_AG_041   Dysbiotic        m_28 24.25   0.000000000     0.00  NA
## MSQ4_AG_042 Normobiotic        m_12 22.84   0.000000000     0.00  NA
## MSQ4_AG_043 Normobiotic        m_13 22.04   0.000000000     0.00  NA
## MSQ4_AG_044 Normobiotic        m_14 22.40   0.000000000     0.00  NA
## MSQ4_AG_045 Normobiotic        m_15 22.25   0.000000000     0.00  NA
## MSQ4_AG_046 Normobiotic        m_25 21.84   0.000000000     0.00  NA
## MSQ4_AG_047         H2O         m_6 21.37   1.340000000    -0.25  NA
## MSQ4_AG_049         H2O        m_31 25.18  -0.710000000    -0.01  NA
## MSQ4_AG_050         H2O        m_33 28.01   2.140000000     0.12  NA
## MSQ4_AG_051         H2O        m_34 26.81  -1.770000000    -0.10  NA
## MSQ4_AG_052      H2O_co        m_10 21.57  -0.240000000    -0.38  NA
## MSQ4_AG_053      H2O_co        m_26 22.52  -2.330000000     0.08  NA
## MSQ4_AG_054      H2O_co        m_27 23.44   2.210000000     0.16  NA
## MSQ4_AG_055      H2O_co        m_29 22.37   0.220000000     0.27  NA
## MSQ4_AG_056      H2O_co        m_30 22.00   1.950000000     0.43  NA
## MSQ4_AG_057         DSS         m_4 20.35  -1.910000000    -0.36  NA
## MSQ4_AG_058         DSS        m_16 22.00  -2.910000000    -0.64  NA
## MSQ4_AG_059         DSS        m_17 22.97   2.480000000     0.59  NA
## MSQ4_AG_062      DSS_co        m_11 23.30  -4.300000000    -0.34  NA
## MSQ4_AG_064      DSS_co        m_22 22.71   1.070000000     0.79  NA
## MSQ4_AG_065      DSS_co        m_23 23.55  -5.330000000     0.68  NA
## MSQ4_AG_072   Dysbiotic         m_7 21.68   0.990000000     0.16  NA
## MSQ4_AG_076   Dysbiotic        m_28 24.98  -0.820000000     0.07  NA
## MSQ4_AG_077 Normobiotic        m_12 21.96  -8.450000000    -0.36  NA
## MSQ4_AG_080 Normobiotic        m_15 22.03  -0.130000000    -0.51  NA
## MSQ4_AG_083         H2O        m_19 23.87   1.420000000    -0.76  NA
## MSQ4_AG_084         H2O        m_31 25.11  -0.080000000    -0.35  NA
## MSQ4_AG_085         H2O        m_33 27.10  -1.960000000    -0.43  NA
## MSQ4_AG_086         H2O        m_34 27.02  -2.800000000    -0.58  NA
## MSQ4_AG_087      H2O_co        m_10 21.96   1.460000000     0.90  NA
## MSQ4_AG_088      H2O_co        m_26 22.74  -0.430000000    -0.47  NA
## MSQ4_AG_089      H2O_co        m_27 23.36  -0.910000000     0.04  NA
## MSQ4_AG_091      H2O_co        m_30 22.24   0.270000000     0.44  NA
## MSQ4_AG_092         DSS        m_16 19.21  -9.450000000    -0.44  NA
## MSQ4_AG_093         DSS        m_17 19.57 -12.790000000    -0.42  NA
## MSQ4_AG_095         DSS        m_20 17.87 -22.910000000    -0.31  NA
## MSQ4_AG_096      DSS_co        m_11 22.55  -2.250000000     0.79  NA
## MSQ4_AG_101    Eubiotic         m_1 22.36  -0.810000000     0.42  NA
## MSQ4_AG_103    Eubiotic         m_3 22.44   0.220000000     0.10  NA
## MSQ4_AG_104    Eubiotic         m_5 22.18   2.450000000     0.09  NA
## MSQ4_AG_105    Eubiotic        m_35 23.18   2.260000000     0.38  NA
## MSQ4_AG_106   Dysbiotic         m_7 19.44 -10.030000000     0.17  NA
## MSQ4_AG_107   Dysbiotic         m_8 19.72 -16.420000000     0.22  NA
## MSQ4_AG_108   Dysbiotic         m_9 21.01  -8.490000000     0.33  NA
## MSQ4_AG_109   Dysbiotic        m_24 20.10  -8.460000000    -0.58  NA
## MSQ4_AG_111 Normobiotic        m_12 20.06 -10.810000000     0.31  NA
## MSQ4_AG_112 Normobiotic        m_13 21.39  -2.590000000     1.27  NA
## MSQ4_AG_113 Normobiotic        m_14 20.84  -7.900000000     1.51  NA
## MSQ4_AG_117         H2O        m_19 24.87   3.670000000     0.32  NA
## MSQ4_AG_119         H2O        m_33 27.12  -1.530000000    -0.73  NA
## MSQ4_AG_121      H2O_co        m_10 22.27   4.850000000     0.58  NA
## MSQ4_AG_125      H2O_co        m_30 22.80   3.640000000     0.12  NA
## MSQ4_AG_127         DSS        m_17 21.20  -7.790000000     0.98  NA
## MSQ4_AG_128         DSS        m_18 20.31 -13.790000000     0.18  NA
## MSQ4_AG_129         DSS        m_20 20.02 -14.440000000     0.35  NA
## MSQ4_AG_130      DSS_co        m_11 24.01   0.210000000    -0.07  NA
## MSQ4_AG_131      DSS_co        m_21 21.05   4.100000000    -0.44  NA
## MSQ4_AG_132      DSS_co        m_22 22.43  -0.180000000    -0.13  NA
## MSQ4_AG_133      DSS_co        m_23 22.16  -5.580000000    -0.55  NA
## MSQ4_AG_135    Eubiotic         m_1 22.58   1.990000000    -0.86  NA
## MSQ4_AG_136    Eubiotic         m_2 22.76   1.380000000    -0.01  NA
## MSQ4_AG_141   Dysbiotic         m_8 20.76 -11.240000000     0.66  NA
## MSQ4_AG_143   Dysbiotic        m_24 20.45  -8.910000000    -0.46  NA
## MSQ4_AG_144   Dysbiotic        m_28 22.73  -6.270000000     0.62  NA
## MSQ4_AG_145 Normobiotic        m_12 20.49 -10.290000000     0.38  NA
## MSQ4_AG_148 Normobiotic        m_15 22.15  -0.450000000     0.27  NA
## MSQ4_AG_149 Normobiotic        m_25 21.93   0.410000000    -0.44  NA
## MSQ4_AG_151         H2O        m_19 24.87   3.670000000     0.32  NA
## MSQ4_AG_152         H2O        m_31 25.42   0.280000000    -0.51  NA
## MSQ4_AG_153         H2O        m_33 27.12  -1.530000000    -0.73  NA
## MSQ4_AG_154         H2O        m_34 27.15  -0.150000000     0.21  NA
## MSQ4_AG_155      H2O_co        m_10 22.27   4.850000000     0.58  NA
## MSQ4_AG_156      H2O_co        m_26 23.10  -0.170000000    -0.08  NA
## MSQ4_AG_157      H2O_co        m_27 23.31   0.950000000     0.13  NA
## MSQ4_AG_159      H2O_co        m_30 22.80   3.640000000     0.12  NA
## MSQ4_AG_160         DSS        m_16 20.65  -6.140000000     0.58  NA
## MSQ4_AG_161         DSS        m_17 21.20  -7.790000000     0.98  NA
## MSQ4_AG_162         DSS        m_18 20.31 -13.790000000     0.18  NA
## MSQ4_AG_163         DSS        m_20 20.02 -14.440000000     0.35  NA
## MSQ4_AG_164      DSS_co        m_11 24.01   0.210000000    -0.07  NA
## MSQ4_AG_165      DSS_co        m_21 21.05   4.100000000    -0.44  NA
## MSQ4_AG_167      DSS_co        m_23 22.16  -5.580000000    -0.55  NA
## MSQ4_AG_168      DSS_co        m_32 22.26  -3.090000000     0.40  NA
## MSQ4_AG_169    Eubiotic         m_1 22.58   1.990000000    -0.86  NA
## MSQ4_AG_170    Eubiotic         m_2 22.76   1.380000000    -0.01  NA
## MSQ4_AG_171    Eubiotic         m_3 22.11  -0.630000000    -0.43  NA
## MSQ4_AG_172    Eubiotic         m_5 23.31   3.920000000     0.53  NA
## MSQ4_AG_173    Eubiotic        m_35 22.82   1.020000000     0.10  NA
## MSQ4_AG_175   Dysbiotic         m_7 20.30  -3.930000000    -0.22  NA
## MSQ4_AG_176   Dysbiotic         m_8 20.76 -11.240000000     0.66  NA
## MSQ4_AG_177   Dysbiotic         m_9 21.85  -4.330000000     0.15  NA
## MSQ4_AG_178   Dysbiotic        m_24 20.45  -8.910000000    -0.46  NA
## MSQ4_AG_179   Dysbiotic        m_28 22.73  -6.270000000     0.62  NA
## MSQ4_AG_180 Normobiotic        m_13 21.82  -1.000000000    -0.16  NA
## MSQ4_AG_181 Normobiotic        m_14 20.80  -7.140000000    -0.06  NA
## MSQ4_AG_183 Normobiotic        m_25 21.93   0.410000000    -0.44  NA
## MSQ4_AG_184         H2O       m_L26 16.70  -6.703910615       NA  NA
## MSQ4_AG_185         H2O       m_L27 16.50  -6.779661017       NA  NA
## MSQ4_AG_186         H2O       m_L28 15.60  -4.878048780       NA  NA
## MSQ4_AG_187         H2O       m_L29 16.30  -9.444444444       NA  NA
## MSQ4_AG_188         H2O       m_L30 12.90  -3.731343284       NA  NA
## MSQ4_AG_189      H2O_co        m_L1 19.30   4.891304348       NA  NA
## MSQ4_AG_191      H2O_co        m_L3 17.80   2.298850575       NA  NA
## MSQ4_AG_192      H2O_co        m_L4 19.10   3.804347826       NA  NA
## MSQ4_AG_193      H2O_co        m_L5 18.70   2.185792350       NA  NA
## MSQ4_AG_194         DSS       m_L31 16.80  -4.000000000       NA  NA
## MSQ4_AG_195         DSS       m_L32 18.10  -1.630434783       NA  NA
## MSQ4_AG_196         DSS       m_L33 16.60  -2.352941176       NA  NA
## MSQ4_AG_197         DSS       m_L34 17.10  -3.932584270       NA  NA
## MSQ4_AG_198         DSS       m_L35 16.40  -2.380952381       NA  NA
## MSQ4_AG_201      DSS_co        m_L9 17.40   1.754385965       NA  NA
## MSQ4_AG_205    Eubiotic       m_L24 18.50  -3.141361257       NA  NA
## MSQ4_AG_207   Dysbiotic       m_L18 17.60   1.149425287       NA  NA
## MSQ4_AG_208   Dysbiotic       m_L19 17.50   0.000000000       NA  NA
## MSQ4_AG_209   Dysbiotic       m_L20 16.80  -1.754385965       NA  NA
## MSQ4_AG_210 Normobiotic       m_L11 18.60  -3.125000000       NA  NA
## MSQ4_AG_211 Normobiotic       m_L12 17.10   1.785714286       NA  NA
## MSQ4_AG_212 Normobiotic       m_L13 19.10  -1.546391753       NA  NA
## MSQ4_AG_213 Normobiotic       m_L15 17.90  -4.787234043       NA  NA
## MSQ4_AG_215         H2O       m_L29 16.30  -9.444444444       NA  NA
## MSQ4_AG_216         H2O       m_L30 12.90  -3.731343284       NA  NA
## MSQ4_AG_217      H2O_co        m_L1 19.30   4.891304348       NA  NA
## MSQ4_AG_218      H2O_co        m_L2 19.60   0.512820513       NA  NA
## MSQ4_AG_219      H2O_co        m_L3 17.80   2.298850575       NA  NA
## MSQ4_AG_220      H2O_co        m_L4 19.10   3.804347826       NA  NA
## MSQ4_AG_221      H2O_co        m_L5 18.70   2.185792350       NA  NA
## MSQ4_AG_223         DSS       m_L32 18.10  -1.630434783       NA  NA
## MSQ4_AG_224         DSS       m_L34 17.10  -3.932584270       NA  NA
## MSQ4_AG_225         DSS       m_L35 16.40  -2.380952381       NA  NA
## MSQ4_AG_226      DSS_co        m_L7 17.10  -3.932584270       NA  NA
## MSQ4_AG_227      DSS_co        m_L9 17.40   1.754385965       NA  NA
## MSQ4_AG_228      DSS_co       m_L10 17.50  -2.234636872       NA  NA
## MSQ4_AG_229    Eubiotic       m_L21 17.90  -1.104972376       NA  NA
## MSQ4_AG_230    Eubiotic       m_L22 16.70   0.602409639       NA  NA
## MSQ4_AG_231    Eubiotic       m_L24 18.50  -3.141361257       NA  NA
## MSQ4_AG_232    Eubiotic       m_L25 17.30  -3.888888889       NA  NA
## MSQ4_AG_233   Dysbiotic       m_L16 17.10  -4.469273743       NA  NA
## MSQ4_AG_234   Dysbiotic       m_L17 16.30  -4.117647059       NA  NA
## MSQ4_AG_235   Dysbiotic       m_L18 17.60   1.149425287       NA  NA
## MSQ4_AG_236   Dysbiotic       m_L19 17.50   0.000000000       NA  NA
## MSQ4_AG_237 Normobiotic       m_L11 18.60  -3.125000000       NA  NA
## MSQ4_AG_238 Normobiotic       m_L12 17.10   1.785714286       NA  NA
## MSQ4_AG_239 Normobiotic       m_L13 19.10  -1.546391753       NA  NA
## MSQ4_AG_240 Normobiotic       m_L15 17.90  -4.787234043       NA  NA
## MSQ4_AG_241         H2O       m_L26 17.90   0.000000000       NA  NA
## MSQ4_AG_242         H2O       m_L27 17.70   0.000000000       NA  NA
## MSQ4_AG_243      H2O_co        m_L1 18.40   0.000000000       NA  NA
## MSQ4_AG_244         DSS       m_L31 17.50   0.000000000       NA  NA
## MSQ4_AG_245      DSS_co        m_L6 17.80   0.000000000       NA  NA
## MSQ4_AG_247   Dysbiotic       m_L16 17.90   0.000000000       NA  NA
## MSQ4_AG_248 Normobiotic       m_L11 19.20   0.000000000       NA  NA
## MSQ4_AG_249 Normobiotic       m_L14 17.90   0.000000000       NA  NA
## MSQ6_AG_250    Eubiotic         m_8 20.55   0.000000000     0.00   0
## MSQ6_AG_251    Eubiotic        m_19 21.32   0.000000000     0.00   0
## MSQ6_AG_252    Eubiotic        m_34 20.11   0.000000000     0.00   0
## MSQ6_AG_253    Eubiotic        m_26 21.90   0.000000000     0.00   0
## MSQ6_AG_254    Eubiotic        m_39 20.94   0.000000000     0.00   0
## MSQ6_AG_255    Eubiotic        m_46 21.12   0.000000000     0.00   0
## MSQ6_AG_256    Eubiotic        m_48 20.35   0.000000000     0.00   0
## MSQ6_AG_257    Eubiotic        m_56 20.51   0.000000000     0.00   0
## MSQ6_AG_258   Dysbiotic        m_29 20.95   0.000000000     0.00   0
## MSQ6_AG_259   Dysbiotic        m_13 20.35   0.000000000     0.00   0
## MSQ6_AG_260   Dysbiotic        m_18 21.55   0.000000000     0.00   0
## MSQ6_AG_261   Dysbiotic        m_27 23.18   0.000000000     0.00   0
## MSQ6_AG_262   Dysbiotic        m_32 20.78   0.000000000     0.00   0
## MSQ6_AG_263   Dysbiotic        m_49 20.90   0.000000000     0.00   0
## MSQ6_AG_264   Dysbiotic        m_54 19.70   0.000000000     0.00   0
## MSQ6_AG_265   Dysbiotic        m_35 20.53   0.000000000     0.00   0
## MSQ6_AG_266 Normobiotic         m_5 22.47   0.000000000     0.00   0
## MSQ6_AG_267 Normobiotic         m_6 21.42   0.000000000     0.00   0
## MSQ6_AG_268 Normobiotic        m_12 20.86   0.000000000     0.00   0
## MSQ6_AG_269 Normobiotic        m_17 21.93   0.000000000     0.00   0
## MSQ6_AG_270 Normobiotic        m_20 21.67   0.000000000     0.00   0
## MSQ6_AG_271 Normobiotic        m_30 21.31   0.000000000     0.00   0
## MSQ6_AG_272 Normobiotic        m_38 21.63   0.000000000     0.00   0
## MSQ6_AG_273 Normobiotic        m_42 21.15   0.000000000     0.00   0
## MSQ6_AG_274      Predni         m_1 20.69   0.000000000     0.00   0
## MSQ6_AG_275      Predni         m_9 21.80   0.000000000     0.00   0
## MSQ6_AG_276      Predni        m_15 22.12   0.000000000     0.00   0
## MSQ6_AG_277      Predni        m_16 20.93   0.000000000     0.00   0
## MSQ6_AG_278      Predni        m_31 22.23   0.000000000     0.00   0
## MSQ6_AG_279      Predni        m_33 21.07   0.000000000     0.00   0
## MSQ6_AG_281      Predni        m_47 21.04   0.000000000     0.00   0
## MSQ6_AG_282         DSS         m_7 20.98   0.000000000     0.00   0
## MSQ6_AG_283         DSS        m_10 20.48   0.000000000     0.00   0
## MSQ6_AG_284         DSS        m_11 20.47   0.000000000     0.00   0
## MSQ6_AG_285         DSS        m_24 22.91   0.000000000     0.00   0
## MSQ6_AG_286         DSS        m_28 21.30   0.000000000     0.00   0
## MSQ6_AG_287         DSS        m_36 22.47   0.000000000     0.00   0
## MSQ6_AG_288         DSS        m_51 21.92   0.000000000     0.00   0
## MSQ6_AG_290         H2O         m_2 21.09   0.000000000     0.00   0
## MSQ6_AG_291         H2O        m_21 21.31   0.000000000     0.00   0
## MSQ6_AG_292         H2O        m_22 20.45   0.000000000     0.00   0
## MSQ6_AG_293         H2O        m_25 22.00   0.000000000     0.00   0
## MSQ6_AG_294         H2O        m_37 20.38   0.000000000     0.00   0
## MSQ6_AG_295         H2O        m_44 21.60   0.000000000     0.00   0
## MSQ6_AG_296         H2O        m_52 21.75   0.000000000     0.00   0
## MSQ6_AG_298    Eubiotic         m_8 18.24   0.112408759     0.20   1
## MSQ6_AG_299    Eubiotic        m_19 18.64   0.125703565     0.64   1
## MSQ6_AG_300    Eubiotic        m_34 19.75   0.017901542     0.90   0
## MSQ6_AG_301    Eubiotic        m_26 22.45  -0.025114155    -0.15   1
## MSQ6_AG_302    Eubiotic        m_39 20.86   0.003820439     1.81   0
## MSQ6_AG_303    Eubiotic        m_46 20.12   0.047348485     0.70   0
## MSQ6_AG_304    Eubiotic        m_48 19.75   0.029484029     0.71   0
## MSQ6_AG_305   Dysbiotic        m_29 17.62   0.158949881     0.26   3
## MSQ6_AG_306   Dysbiotic        m_13 17.60   0.135135135     0.05   3
## MSQ6_AG_307   Dysbiotic        m_18 17.57   0.184686775     0.39   3
## MSQ6_AG_308   Dysbiotic        m_27 20.47   0.116911130    -0.46   1
## MSQ6_AG_309   Dysbiotic        m_32 17.32   0.166506256     0.24   3
## MSQ6_AG_310   Dysbiotic        m_49 17.95   0.141148325     0.18   3
## MSQ6_AG_311   Dysbiotic        m_54 17.25   0.124365482     0.20   3
## MSQ6_AG_312   Dysbiotic        m_35 18.45   0.101315149    -0.10   3
## MSQ6_AG_313 Normobiotic         m_5 18.68   0.168669337     0.62   2
## MSQ6_AG_314 Normobiotic         m_6 18.90   0.117647059     0.52   2
## MSQ6_AG_315 Normobiotic        m_12 19.36   0.071907958     0.76   1
## MSQ6_AG_316 Normobiotic        m_17 20.69   0.056543548     1.14   1
## MSQ6_AG_317 Normobiotic        m_20 18.89   0.128287956     0.86   3
## MSQ6_AG_318 Normobiotic        m_30 18.00   0.155326138     0.15   3
## MSQ6_AG_319 Normobiotic        m_38 19.59   0.094313454     1.08   1
## MSQ6_AG_320 Normobiotic        m_42 19.72   0.067612293     0.29   2
## MSQ6_AG_321      Predni         m_1 17.51   0.153697438     1.05   4
## MSQ6_AG_322      Predni         m_9 20.32   0.067889908     1.20   1
## MSQ6_AG_323      Predni        m_15 20.37   0.079113924     0.96   2
## MSQ6_AG_324      Predni        m_16 17.39   0.169135213     0.93   3
## MSQ6_AG_325      Predni        m_31 19.27   0.133153396     0.20   3
## MSQ6_AG_326      Predni        m_33 16.81   0.202183199     0.41   4
## MSQ6_AG_327      Predni        m_40 18.14   0.180298238     1.06   4
## MSQ6_AG_328      Predni        m_47 18.48   0.121673004     1.15   2
## MSQ6_AG_329         DSS         m_7 18.25   0.130123928     0.30   3
## MSQ6_AG_330         DSS        m_10 17.35   0.152832031     0.25   3
## MSQ6_AG_331         DSS        m_11 16.85   0.176844162    -0.01   4
## MSQ6_AG_332         DSS        m_24 20.01   0.126582278    -0.35   2
## MSQ6_AG_333         DSS        m_28 18.12   0.149295775    -0.05   3
## MSQ6_AG_334         DSS        m_36 17.07   0.240320427    -0.19   4
## MSQ6_AG_335         DSS        m_51 16.07   0.266879562     0.20   4
## MSQ6_AG_336         H2O         m_2 20.70  -0.018492176    -0.02   0
## MSQ6_AG_337         H2O        m_21 21.54   0.010793055     0.74   0
## MSQ6_AG_339         H2O        m_25 21.86  -0.006363636    -1.88   0
## MSQ6_AG_340         H2O        m_37 21.16   0.038272816    -0.22   0
## MSQ6_AG_341         H2O        m_44 22.97   0.063425926     0.12   0
## MSQ6_AG_342         H2O        m_52 22.07   0.014712644     0.67   0
## MSQ6_AG_343         H2O        m_53 21.47  -0.029823769    -0.87   0
## MSQ6_AG_344    Eubiotic        m_19 21.13   0.008911820     1.32   0
## MSQ6_AG_346    Eubiotic        m_34 21.03  -0.045748384     0.84   0
## MSQ6_AG_347    Eubiotic        m_26 21.34   0.025570776    -0.85   0
## MSQ6_AG_348    Eubiotic        m_39 22.88  -0.092645654     2.42   0
## MSQ6_AG_349    Eubiotic        m_46 21.81  -0.032670455     0.59   0
## MSQ6_AG_350    Eubiotic        m_48 21.31  -0.047174447     0.49   0
## MSQ6_AG_351   Dysbiotic        m_29 18.02   0.139856802     0.17   2
## MSQ6_AG_352   Dysbiotic        m_13 18.23   0.104176904     0.38   2
## MSQ6_AG_353   Dysbiotic        m_18 18.01   0.164269142     0.45   3
## MSQ6_AG_354   Dysbiotic        m_27 20.95   0.096203624     0.94   0
## MSQ6_AG_355   Dysbiotic        m_32 17.85   0.141000962     0.50   2
## MSQ6_AG_356   Dysbiotic        m_49 18.35   0.122009569     0.23   2
## MSQ6_AG_357   Dysbiotic        m_54 18.52   0.059898477     0.52   0
## MSQ6_AG_358   Dysbiotic        m_35 18.95   0.076960546     0.41   0
## MSQ6_AG_359 Normobiotic         m_5 20.02   0.109034268     0.34   2
## MSQ6_AG_360 Normobiotic         m_6 19.65   0.082633053     0.65   1
## MSQ6_AG_361 Normobiotic        m_12 20.53   0.015819751     1.00   0
## MSQ6_AG_362 Normobiotic        m_17 21.02   0.041495668     0.33   0
## MSQ6_AG_363 Normobiotic        m_20 19.93   0.080295339     0.40   0
## MSQ6_AG_364 Normobiotic        m_42 20.23   0.043498818     0.43   1
## MSQ6_AG_365      Predni         m_1 18.41   0.110198163    -0.15   2
## MSQ6_AG_366      Predni         m_9 21.20   0.027522936    -0.17   1
## MSQ6_AG_367      Predni        m_15 21.76   0.016274864     0.47   0
## MSQ6_AG_368      Predni        m_31 20.50   0.077822762     0.43   1
## MSQ6_AG_369      Predni        m_33 17.85   0.152823920     1.33   4
## MSQ6_AG_370      Predni        m_40 20.70   0.064618165     1.08   1
## MSQ6_AG_371      Predni        m_47 20.20   0.039923954     0.97   0
## MSQ6_AG_372         DSS         m_7 18.21   0.132030505     0.07   2
## MSQ6_AG_373         DSS        m_10 17.58   0.141601563    -0.07   3
## MSQ6_AG_374         DSS        m_11 17.65   0.137762579     0.09   2
## MSQ6_AG_375         DSS        m_24 20.23   0.116979485    -0.40   3
## MSQ6_AG_376         DSS        m_28 18.96   0.109859155     0.42   3
## MSQ6_AG_377         DSS        m_36 18.02   0.198041834     0.80   3
## MSQ6_AG_378         DSS        m_51 16.96   0.226277372     0.72   4
## MSQ6_AG_379         DSS        m_55 18.02   0.121404193     0.06   2
## MSQ6_AG_380         H2O         m_2 21.28   0.009009009     0.58   0
## MSQ6_AG_381         H2O        m_21 21.72   0.019239794     0.18   0
## MSQ6_AG_382         H2O        m_22 21.70   0.061124694     0.54   0
## MSQ6_AG_383         H2O        m_25 21.79  -0.009545455    -0.07   0
## MSQ6_AG_384         H2O        m_37 21.72   0.065750736     0.56   0
## MSQ6_AG_385         H2O        m_44 23.77   0.100462963     0.80   0
## MSQ6_AG_386         H2O        m_52 21.75   0.000000000    -0.26   0
## MSQ6_AG_387         H2O        m_53 22.26   0.005874379     0.79   0
## MSQ6_AG_388    Eubiotic         m_8 20.21   0.016545012     0.63   0
## MSQ6_AG_389    Eubiotic        m_19 21.13   0.008911820     1.32   0
## MSQ6_AG_390    Eubiotic        m_34 21.03  -0.045748384     0.84   0
## MSQ6_AG_391    Eubiotic        m_26 21.34   0.025570776    -0.85   0
## MSQ6_AG_392    Eubiotic        m_39 22.88  -0.092645654     2.42   0
## MSQ6_AG_393    Eubiotic        m_46 21.81  -0.032670455     0.59   0
## MSQ6_AG_394    Eubiotic        m_48 21.31  -0.047174447     0.49   0
## MSQ6_AG_395   Dysbiotic        m_29 18.02   0.139856802     0.17   2
## MSQ6_AG_396   Dysbiotic        m_13 18.23   0.104176904     0.38   2
## MSQ6_AG_397   Dysbiotic        m_18 18.01   0.164269142     0.45   3
## MSQ6_AG_398   Dysbiotic        m_27 20.95   0.096203624     0.94   0
## MSQ6_AG_399   Dysbiotic        m_32 17.85   0.141000962     0.50   2
## MSQ6_AG_400   Dysbiotic        m_49 18.35   0.122009569     0.23   0
## MSQ6_AG_401   Dysbiotic        m_54 18.52   0.059898477     0.52   2
## MSQ6_AG_402   Dysbiotic        m_35 18.95   0.076960546     0.41   0
## MSQ6_AG_403 Normobiotic         m_5 20.02   0.109034268     0.34   1
## MSQ6_AG_404 Normobiotic         m_6 19.65   0.082633053     0.65   0
## MSQ6_AG_405 Normobiotic        m_12 20.53   0.015819751     1.00   0
## MSQ6_AG_406 Normobiotic        m_17 21.02   0.041495668     0.33   0
## MSQ6_AG_407 Normobiotic        m_20 19.93   0.080295339     0.40   0
## MSQ6_AG_408 Normobiotic        m_30 19.35   0.091975598     0.46   1
## MSQ6_AG_409 Normobiotic        m_38 20.56   0.049468331     0.43   1
## MSQ6_AG_410 Normobiotic        m_42 20.23   0.043498818     0.43   1
## MSQ6_AG_411      Predni         m_1 18.41   0.110198163    -0.15   2
## MSQ6_AG_412      Predni         m_9 21.20   0.027522936    -0.17   1
## MSQ6_AG_413      Predni        m_15 21.76   0.016274864     0.47   0
## MSQ6_AG_414      Predni        m_16 18.03   0.138557095     0.42   3
## MSQ6_AG_415      Predni        m_31 20.50   0.077822762     0.43   1
## MSQ6_AG_416      Predni        m_33 17.85   0.152823920     1.33   4
## MSQ6_AG_417      Predni        m_40 20.70   0.064618165     1.08   1
## MSQ6_AG_418      Predni        m_47 20.20   0.039923954     0.97   0
## MSQ6_AG_419         DSS         m_7 18.21   0.132030505     0.07   2
## MSQ6_AG_420         DSS        m_10 17.58   0.141601563    -0.07   3
## MSQ6_AG_421         DSS        m_11 17.65   0.137762579     0.09   2
## MSQ6_AG_422         DSS        m_24 20.23   0.116979485    -0.40   3
## MSQ6_AG_423         DSS        m_28 18.96   0.109859155     0.42   3
## MSQ6_AG_424         DSS        m_36 18.02   0.198041834     0.80   3
## MSQ6_AG_425         DSS        m_51 16.96   0.226277372     0.72   4
## MSQ6_AG_426         DSS        m_55 18.02   0.121404193     0.06   2
## MSQ6_AG_427         H2O         m_2 21.28   0.009009009     0.58   0
## MSQ6_AG_428         H2O        m_21 21.72   0.019239794     0.18   0
## MSQ6_AG_429         H2O        m_22 21.70   0.061124694     0.54   0
## MSQ6_AG_430         H2O        m_25 21.79  -0.009545455    -0.07   0
## MSQ6_AG_431         H2O        m_37 21.72   0.065750736     0.56   0
## MSQ6_AG_432         H2O        m_44 23.77   0.100462963     0.80   0
## MSQ6_AG_433         H2O        m_52 21.75   0.000000000    -0.26   0
## MSQ6_AG_434         H2O        m_53 22.26   0.005874379     0.79   0
## MSQ6_AG_490    Eubiotic        m_19 20.48   0.039399625    -0.93   1
## MSQ6_AG_491    Eubiotic        m_26 20.35   0.070776256    -1.61   2
## MSQ6_AG_492    Eubiotic        m_46 20.77   0.016571970    -0.49   1
## MSQ6_AG_493    Eubiotic        m_48 19.90   0.022113022    -0.43   1
## MSQ6_AG_494    Eubiotic        m_56 19.50   0.049244271    -1.12   1
## MSQ6_AG_495   Dysbiotic        m_29 21.32  -0.017661098     0.39   2
## MSQ6_AG_496   Dysbiotic        m_13 20.59  -0.011793612    -0.79   1
## MSQ6_AG_497   Dysbiotic        m_18 21.06   0.022737819    -1.04   1
## MSQ6_AG_498   Dysbiotic        m_27 22.95   0.009922347    -0.21   0
## MSQ6_AG_499   Dysbiotic        m_32 20.91  -0.006256015    -0.15   1
## MSQ6_AG_500   Dysbiotic        m_49 20.41   0.023444976    -1.00   2
## MSQ6_AG_501   Dysbiotic        m_54 19.93  -0.011675127    -0.86   1
## MSQ6_AG_502   Dysbiotic        m_35 20.81  -0.013638578    -0.60   1
## MSQ6_AG_503 Normobiotic         m_5 20.76   0.076101469    -1.09   1
## MSQ6_AG_504 Normobiotic         m_6 21.88  -0.021475257     0.19   1
## MSQ6_AG_505 Normobiotic        m_12 20.48   0.018216683    -0.37   1
## MSQ6_AG_506 Normobiotic        m_17 22.09  -0.007295942    -0.31   0
## MSQ6_AG_509 Normobiotic        m_38 21.92  -0.013407305    -0.46   0
## MSQ6_AG_511      Predni        m_15 22.13  -0.000452000     0.58   1
## MSQ6_AG_512      Predni        m_16 20.07   0.041089345    -0.36   1
## MSQ6_AG_514      Predni        m_33 19.75   0.062648315    -0.91   2
## MSQ6_AG_515      Predni        m_40 21.20   0.042024401    -0.78   2
## MSQ6_AG_516      Predni        m_47 20.80   0.011406844    -0.09   1
## MSQ6_AG_517         DSS         m_7 20.98   0.000000000    -0.93   1
## MSQ6_AG_518         DSS        m_10 20.42   0.002929687    -0.69   1
## MSQ6_AG_519         DSS        m_11 19.93   0.026380068    -0.87   1
## MSQ6_AG_520         DSS        m_24 23.81  -0.039284155     0.46   0
## MSQ6_AG_521         DSS        m_28 20.73   0.026760563    -0.70   1
## MSQ6_AG_522         DSS        m_36 21.79   0.030262572    -0.33   2
## MSQ6_AG_523         H2O         m_2 20.65  -0.020862968    -0.34   0
## MSQ6_AG_524         H2O        m_21 21.20  -0.005161896     0.25   0
## MSQ6_AG_525         H2O        m_22 20.90   0.022004890    -1.05   0
## MSQ6_AG_526         H2O        m_25 21.68  -0.014545455     0.38   0
## MSQ6_AG_527         H2O        m_37 21.68   0.063788027     0.43   0
## MSQ6_AG_528         H2O        m_44 22.45   0.039351852     0.52   0
## MSQ6_AG_529         H2O        m_52 20.76  -0.045517241    -0.19   0
## MSQ6_AG_530         H2O        m_53 21.74  -0.017623136     0.18   0
##             colon_lenght spleen_weight histo_colo_infilt histo_colo_damage
## MSQ4_AG_012          5.5        0.1510               0.5              0.50
## MSQ4_AG_013          5.4        0.0860               0.0              0.00
## MSQ4_AG_014          5.7        0.1300               0.0              0.00
## MSQ4_AG_015          5.8        0.1160               0.0              0.00
## MSQ4_AG_016          6.3        0.1330               0.5              0.50
## MSQ4_AG_017          5.6        0.1210               0.0              0.00
## MSQ4_AG_018          5.7        0.1600               0.0              0.00
## MSQ4_AG_019          5.4        0.1380               0.5              0.50
## MSQ4_AG_020          6.1        0.1270               0.0              0.00
## MSQ4_AG_021          6.0        0.1240               0.0              0.00
## MSQ4_AG_022           NA            NA                NA                NA
## MSQ4_AG_023          5.3        0.1420               4.0              4.00
## MSQ4_AG_024          5.7        0.1800               4.0              3.50
## MSQ4_AG_025          5.7        0.1820               3.5              3.50
## MSQ4_AG_027          6.4        0.1200               1.5              2.00
## MSQ4_AG_028          7.4        0.0920               1.5              2.00
## MSQ4_AG_029          7.3        0.1150               1.5              1.50
## MSQ4_AG_030          6.8        0.1280               1.0              1.50
## MSQ4_AG_031          7.2        0.0890               1.0              1.00
## MSQ4_AG_033          8.2        0.0950               3.0              3.00
## MSQ4_AG_034          8.3        0.0940               0.0              0.00
## MSQ4_AG_035          9.1        0.0980               2.5              2.00
## MSQ4_AG_036          7.8        0.1020               1.5              1.00
## MSQ4_AG_037          8.2        0.0840               4.0              3.50
## MSQ4_AG_038          7.9        0.0760               4.0              4.00
## MSQ4_AG_039          8.2        0.0910               3.5              3.00
## MSQ4_AG_040          7.8        0.0820               3.0              3.00
## MSQ4_AG_041          7.9        0.0860               4.0              4.00
## MSQ4_AG_042          6.9        0.0910               4.0              3.75
## MSQ4_AG_043          5.8        0.1010               0.0              0.50
## MSQ4_AG_044          5.9        0.1070               4.0              3.50
## MSQ4_AG_045          6.2        0.0980               1.0              1.00
## MSQ4_AG_046          6.8        0.1120               3.0              3.00
## MSQ4_AG_047          5.5        0.1510               0.5              0.50
## MSQ4_AG_049          5.7        0.1300               0.0              0.00
## MSQ4_AG_050          5.8        0.1160               0.0              0.00
## MSQ4_AG_051          6.3        0.1330               0.5              0.50
## MSQ4_AG_052          5.6        0.1210               0.0              0.00
## MSQ4_AG_053          5.7        0.1600               0.0              0.00
## MSQ4_AG_054          5.4        0.1380               0.5              0.50
## MSQ4_AG_055          6.1        0.1270               0.0              0.00
## MSQ4_AG_056          6.0        0.1240               0.0              0.00
## MSQ4_AG_057           NA            NA                NA                NA
## MSQ4_AG_058          5.3        0.1420               4.0              4.00
## MSQ4_AG_059          5.7        0.1800               4.0              3.50
## MSQ4_AG_062          6.4        0.1200               1.5              2.00
## MSQ4_AG_064          7.3        0.1150               1.5              1.50
## MSQ4_AG_065          6.8        0.1280               1.0              1.50
## MSQ4_AG_072          8.2        0.0840               4.0              3.50
## MSQ4_AG_076          7.9        0.0860               4.0              4.00
## MSQ4_AG_077          6.9        0.0910               4.0              3.75
## MSQ4_AG_080          6.2        0.1010               1.0              1.00
## MSQ4_AG_083          5.4        0.0860               0.0              0.00
## MSQ4_AG_084          5.7        0.1300               0.0              0.00
## MSQ4_AG_085          5.8        0.1160               0.0              0.00
## MSQ4_AG_086          6.3        0.1330               0.5              0.50
## MSQ4_AG_087          5.6        0.1210               0.0              0.00
## MSQ4_AG_088          5.7        0.1600               0.0              0.00
## MSQ4_AG_089          5.4        0.1380               0.5              0.50
## MSQ4_AG_091          6.0        0.1240               0.0              0.00
## MSQ4_AG_092          5.3        0.1420               4.0              4.00
## MSQ4_AG_093          5.7        0.1800               4.0              3.50
## MSQ4_AG_095          5.9        0.1790               3.0              3.50
## MSQ4_AG_096          6.4        0.1200               1.5              2.00
## MSQ4_AG_101          7.5        0.0850               1.5              2.50
## MSQ4_AG_103          8.3        0.0940               0.0              0.00
## MSQ4_AG_104          9.1        0.0980               2.5              2.00
## MSQ4_AG_105          7.8        0.1020               1.5              1.00
## MSQ4_AG_106          8.2        0.0840               4.0              3.50
## MSQ4_AG_107          7.9        0.0760               4.0              4.00
## MSQ4_AG_108          8.2        0.0910               3.5              3.00
## MSQ4_AG_109          7.8        0.0820               3.0              3.00
## MSQ4_AG_111          6.9        0.0910               4.0              3.75
## MSQ4_AG_112          5.8        0.1010               0.0              0.50
## MSQ4_AG_113          5.9        0.1070               4.0              3.50
## MSQ4_AG_117          5.4        0.0860               0.0              0.00
## MSQ4_AG_119          5.8        0.1160               0.0              0.00
## MSQ4_AG_121          5.6        0.1210               0.0              0.00
## MSQ4_AG_125          6.0        0.1240               0.0              0.00
## MSQ4_AG_127          5.7        0.1800               4.0              3.50
## MSQ4_AG_128          5.7        0.1820               3.5              3.50
## MSQ4_AG_129          5.9        0.1790               3.0              3.50
## MSQ4_AG_130          6.4        0.1200               1.5              2.00
## MSQ4_AG_131          7.4        0.0920               1.5              2.00
## MSQ4_AG_132          7.3        0.1150               1.5              1.50
## MSQ4_AG_133          6.8        0.1280               1.0              1.50
## MSQ4_AG_135          7.5        0.0850               1.5              2.50
## MSQ4_AG_136          8.2        0.0950               3.0              3.00
## MSQ4_AG_141          7.9        0.0760               4.0              4.00
## MSQ4_AG_143          7.8        0.0820               3.0              3.00
## MSQ4_AG_144          7.9        0.0860               4.0              4.00
## MSQ4_AG_145          6.9        0.0910               4.0              3.75
## MSQ4_AG_148          6.2        0.0980               1.0              1.00
## MSQ4_AG_149          6.8        0.1120               3.0              3.00
## MSQ4_AG_151          5.4        0.0860               0.0              0.00
## MSQ4_AG_152          5.7        0.1300               0.0              0.00
## MSQ4_AG_153          5.8        0.1160               0.0              0.00
## MSQ4_AG_154          6.3        0.1330               0.5              0.50
## MSQ4_AG_155          5.6        0.1210               0.0              0.00
## MSQ4_AG_156          5.7        0.1600               0.0              0.00
## MSQ4_AG_157          5.4        0.1380               0.5              0.50
## MSQ4_AG_159          6.0        0.1240               0.0              0.00
## MSQ4_AG_160          5.3        0.1420               4.0              4.00
## MSQ4_AG_161          5.7        0.1800               4.0              3.50
## MSQ4_AG_162          5.7        0.1820               3.5              3.50
## MSQ4_AG_163          5.9        0.1790               3.0              3.50
## MSQ4_AG_164          6.4        0.1200               1.5              2.00
## MSQ4_AG_165          7.4        0.0920               1.5              2.00
## MSQ4_AG_167          6.8        0.1280               1.0              1.50
## MSQ4_AG_168          7.2        0.0890               1.0              1.00
## MSQ4_AG_169          7.5        0.0850               1.5              2.50
## MSQ4_AG_170          8.2        0.0950               3.0              3.00
## MSQ4_AG_171          8.3        0.0940               0.0              0.00
## MSQ4_AG_172          9.1        0.0980               2.5              2.00
## MSQ4_AG_173          7.8        0.1020               1.5              1.00
## MSQ4_AG_175          8.2        0.0840               4.0              3.50
## MSQ4_AG_176          7.9        0.0760               4.0              4.00
## MSQ4_AG_177          8.2        0.0910               3.5              3.00
## MSQ4_AG_178          7.8        0.0820               3.0              3.00
## MSQ4_AG_179          7.9        0.0860               4.0              4.00
## MSQ4_AG_180          5.8        0.1010               0.0              0.50
## MSQ4_AG_181          5.9        0.1070               4.0              3.50
## MSQ4_AG_183          6.8        0.1120               3.0              3.00
## MSQ4_AG_184          7.1        0.1077                NA                NA
## MSQ4_AG_185          6.9        0.1338                NA                NA
## MSQ4_AG_186          7.6        0.1016                NA                NA
## MSQ4_AG_187          7.4        0.1139                NA                NA
## MSQ4_AG_188          6.5        0.0546                NA                NA
## MSQ4_AG_189          8.3        0.0981                NA                NA
## MSQ4_AG_191          7.7        0.0995                NA                NA
## MSQ4_AG_192          8.1        0.1078                NA                NA
## MSQ4_AG_193          8.2        0.1182                NA                NA
## MSQ4_AG_194          6.3        0.1602                NA                NA
## MSQ4_AG_195          7.1        0.1604                NA                NA
## MSQ4_AG_196          6.7        0.1625                NA                NA
## MSQ4_AG_197          6.2        0.1768                NA                NA
## MSQ4_AG_198          6.1        0.2168                NA                NA
## MSQ4_AG_201          5.3        0.1668                NA                NA
## MSQ4_AG_205          5.1        0.2205                NA                NA
## MSQ4_AG_207          4.1        0.1776                NA                NA
## MSQ4_AG_208          4.9        0.1749                NA                NA
## MSQ4_AG_209          5.5        0.1484                NA                NA
## MSQ4_AG_210          5.2        0.1921                NA                NA
## MSQ4_AG_211          5.4        0.1508                NA                NA
## MSQ4_AG_212          5.4        0.1889                NA                NA
## MSQ4_AG_213          5.1        0.1965                NA                NA
## MSQ4_AG_215          7.4        0.1139                NA                NA
## MSQ4_AG_216          6.5        0.0546                NA                NA
## MSQ4_AG_217          8.3        0.0981                NA                NA
## MSQ4_AG_218          7.4        0.1028                NA                NA
## MSQ4_AG_219          7.7        0.0995                NA                NA
## MSQ4_AG_220          8.1        0.1078                NA                NA
## MSQ4_AG_221          8.2        0.1182                NA                NA
## MSQ4_AG_223          7.1        0.1604                NA                NA
## MSQ4_AG_224          6.2        0.1768                NA                NA
## MSQ4_AG_225          6.1        0.2168                NA                NA
## MSQ4_AG_226          6.0        0.1664                NA                NA
## MSQ4_AG_227          5.3        0.1668                NA                NA
## MSQ4_AG_228          6.1        0.1536                NA                NA
## MSQ4_AG_229          5.8        0.2200                NA                NA
## MSQ4_AG_230          5.6        0.1830                NA                NA
## MSQ4_AG_231          5.1        0.2205                NA                NA
## MSQ4_AG_232          5.4        0.1862                NA                NA
## MSQ4_AG_233          4.9        0.1429                NA                NA
## MSQ4_AG_234          5.5        0.1506                NA                NA
## MSQ4_AG_235          4.1        0.1776                NA                NA
## MSQ4_AG_236          4.9        0.1749                NA                NA
## MSQ4_AG_237          5.2        0.1921                NA                NA
## MSQ4_AG_238          5.4        0.1508                NA                NA
## MSQ4_AG_239          5.4        0.1889                NA                NA
## MSQ4_AG_240          5.1        0.1965                NA                NA
## MSQ4_AG_241          7.1        0.1077                NA                NA
## MSQ4_AG_242          6.9        0.1338                NA                NA
## MSQ4_AG_243          8.3        0.0981                NA                NA
## MSQ4_AG_244          6.3        0.1602                NA                NA
## MSQ4_AG_245          5.0        0.0382                NA                NA
## MSQ4_AG_247          4.9        0.1429                NA                NA
## MSQ4_AG_248          5.2        0.1921                NA                NA
## MSQ4_AG_249          5.0        0.0493                NA                NA
## MSQ6_AG_250           NA            NA                NA                NA
## MSQ6_AG_251           NA            NA                NA                NA
## MSQ6_AG_252           NA            NA                NA                NA
## MSQ6_AG_253           NA            NA                NA                NA
## MSQ6_AG_254           NA            NA                NA                NA
## MSQ6_AG_255           NA            NA                NA                NA
## MSQ6_AG_256           NA            NA                NA                NA
## MSQ6_AG_257           NA            NA                NA                NA
## MSQ6_AG_258           NA            NA                NA                NA
## MSQ6_AG_259           NA            NA                NA                NA
## MSQ6_AG_260           NA            NA                NA                NA
## MSQ6_AG_261           NA            NA                NA                NA
## MSQ6_AG_262           NA            NA                NA                NA
## MSQ6_AG_263           NA            NA                NA                NA
## MSQ6_AG_264           NA            NA                NA                NA
## MSQ6_AG_265           NA            NA                NA                NA
## MSQ6_AG_266           NA            NA                NA                NA
## MSQ6_AG_267           NA            NA                NA                NA
## MSQ6_AG_268           NA            NA                NA                NA
## MSQ6_AG_269           NA            NA                NA                NA
## MSQ6_AG_270           NA            NA                NA                NA
## MSQ6_AG_271           NA            NA                NA                NA
## MSQ6_AG_272           NA            NA                NA                NA
## MSQ6_AG_273           NA            NA                NA                NA
## MSQ6_AG_274           NA            NA                NA                NA
## MSQ6_AG_275           NA            NA                NA                NA
## MSQ6_AG_276           NA            NA                NA                NA
## MSQ6_AG_277           NA            NA                NA                NA
## MSQ6_AG_278           NA            NA                NA                NA
## MSQ6_AG_279           NA            NA                NA                NA
## MSQ6_AG_281           NA            NA                NA                NA
## MSQ6_AG_282           NA            NA                NA                NA
## MSQ6_AG_283           NA            NA                NA                NA
## MSQ6_AG_284           NA            NA                NA                NA
## MSQ6_AG_285           NA            NA                NA                NA
## MSQ6_AG_286           NA            NA                NA                NA
## MSQ6_AG_287           NA            NA                NA                NA
## MSQ6_AG_288           NA            NA                NA                NA
## MSQ6_AG_290           NA            NA                NA                NA
## MSQ6_AG_291           NA            NA                NA                NA
## MSQ6_AG_292           NA            NA                NA                NA
## MSQ6_AG_293           NA            NA                NA                NA
## MSQ6_AG_294           NA            NA                NA                NA
## MSQ6_AG_295           NA            NA                NA                NA
## MSQ6_AG_296           NA            NA                NA                NA
## MSQ6_AG_298           NA            NA                NA                NA
## MSQ6_AG_299           NA            NA                NA                NA
## MSQ6_AG_300           NA            NA                NA                NA
## MSQ6_AG_301           NA            NA                NA                NA
## MSQ6_AG_302           NA            NA                NA                NA
## MSQ6_AG_303           NA            NA                NA                NA
## MSQ6_AG_304           NA            NA                NA                NA
## MSQ6_AG_305           NA            NA                NA                NA
## MSQ6_AG_306           NA            NA                NA                NA
## MSQ6_AG_307           NA            NA                NA                NA
## MSQ6_AG_308           NA            NA                NA                NA
## MSQ6_AG_309           NA            NA                NA                NA
## MSQ6_AG_310           NA            NA                NA                NA
## MSQ6_AG_311           NA            NA                NA                NA
## MSQ6_AG_312           NA            NA                NA                NA
## MSQ6_AG_313           NA            NA                NA                NA
## MSQ6_AG_314           NA            NA                NA                NA
## MSQ6_AG_315           NA            NA                NA                NA
## MSQ6_AG_316           NA            NA                NA                NA
## MSQ6_AG_317           NA            NA                NA                NA
## MSQ6_AG_318           NA            NA                NA                NA
## MSQ6_AG_319           NA            NA                NA                NA
## MSQ6_AG_320           NA            NA                NA                NA
## MSQ6_AG_321           NA            NA                NA                NA
## MSQ6_AG_322           NA            NA                NA                NA
## MSQ6_AG_323           NA            NA                NA                NA
## MSQ6_AG_324           NA            NA                NA                NA
## MSQ6_AG_325           NA            NA                NA                NA
## MSQ6_AG_326           NA            NA                NA                NA
## MSQ6_AG_327           NA            NA                NA                NA
## MSQ6_AG_328           NA            NA                NA                NA
## MSQ6_AG_329           NA            NA                NA                NA
## MSQ6_AG_330           NA            NA                NA                NA
## MSQ6_AG_331           NA            NA                NA                NA
## MSQ6_AG_332           NA            NA                NA                NA
## MSQ6_AG_333           NA            NA                NA                NA
## MSQ6_AG_334           NA            NA                NA                NA
## MSQ6_AG_335           NA            NA                NA                NA
## MSQ6_AG_336           NA            NA                NA                NA
## MSQ6_AG_337           NA            NA                NA                NA
## MSQ6_AG_339           NA            NA                NA                NA
## MSQ6_AG_340           NA            NA                NA                NA
## MSQ6_AG_341           NA            NA                NA                NA
## MSQ6_AG_342           NA            NA                NA                NA
## MSQ6_AG_343           NA            NA                NA                NA
## MSQ6_AG_344          5.2        0.1050               1.0              1.50
## MSQ6_AG_346          5.4        0.0860               2.5              2.00
## MSQ6_AG_347          6.0        0.0950               2.5              2.00
## MSQ6_AG_348          5.2        0.1010               1.0              1.00
## MSQ6_AG_349          6.1        0.1280               1.5              1.00
## MSQ6_AG_350          5.8        0.1010               2.0              2.00
## MSQ6_AG_351          4.0        0.1540               3.0              3.50
## MSQ6_AG_352          4.2        0.1620               4.0              4.00
## MSQ6_AG_353          4.5        0.1440               3.0              2.50
## MSQ6_AG_354          4.2        0.1280               3.5              3.00
## MSQ6_AG_355          4.6        0.1520               2.5              2.50
## MSQ6_AG_356          4.9        0.1290               3.0              4.00
## MSQ6_AG_357          5.1        0.1360               3.0              3.00
## MSQ6_AG_358          4.3        0.1480               1.5              2.00
## MSQ6_AG_359          4.9        0.1670               3.0              2.50
## MSQ6_AG_360          5.3        0.1180               2.5              3.00
## MSQ6_AG_361          5.8        0.1150               2.0              2.50
## MSQ6_AG_362          4.9        0.1240               3.0              3.00
## MSQ6_AG_363          5.2        0.1400               2.5              3.00
## MSQ6_AG_364          5.1        0.1140               1.0              2.00
## MSQ6_AG_365          5.4        0.1030               2.5              3.00
## MSQ6_AG_366          5.2        0.1430               4.0              4.00
## MSQ6_AG_367          4.9        0.0960               3.0              3.50
## MSQ6_AG_368          4.5        0.0960               2.5              3.00
## MSQ6_AG_369          5.1        0.0860               2.0              2.50
## MSQ6_AG_370          4.8        0.1300               1.5              2.00
## MSQ6_AG_371          4.7        0.1200               3.0              4.00
## MSQ6_AG_372          4.2        0.1470               4.0              4.00
## MSQ6_AG_373          4.7        0.1230               3.5              3.50
## MSQ6_AG_374          4.3        0.1250               3.0              4.00
## MSQ6_AG_375          4.6        0.1420               2.5              3.00
## MSQ6_AG_376          4.4        0.1100               4.0              4.00
## MSQ6_AG_377          4.9        0.1830               3.5              4.00
## MSQ6_AG_378          4.3        0.1650               3.5              4.00
## MSQ6_AG_379          4.2        0.1450               4.0              3.00
## MSQ6_AG_380          7.8        0.0650               0.5              0.00
## MSQ6_AG_381          6.2        0.0930               0.0              0.00
## MSQ6_AG_382          6.9        0.0740               1.0              0.50
## MSQ6_AG_383          7.1        0.0710               0.0              0.00
## MSQ6_AG_384          6.8        0.0670               0.0              0.00
## MSQ6_AG_385          6.7        0.0820               0.0              0.00
## MSQ6_AG_386          6.5        0.0650               0.0              0.00
## MSQ6_AG_387          6.9        0.0650               0.0              0.50
## MSQ6_AG_388          5.7        0.1240               1.5              2.00
## MSQ6_AG_389          5.2        0.1050               1.0              1.50
## MSQ6_AG_390          5.4        0.0860               2.5              2.00
## MSQ6_AG_391          6.0        0.0950               2.5              2.00
## MSQ6_AG_392          5.2        0.1010               1.0              1.00
## MSQ6_AG_393          6.1        0.1280               1.5              1.00
## MSQ6_AG_394          5.8        0.1010               2.0              2.00
## MSQ6_AG_395          4.0        0.1540               3.0              3.50
## MSQ6_AG_396          4.2        0.1620               4.0              4.00
## MSQ6_AG_397          4.5        0.1440               3.0              2.50
## MSQ6_AG_398          4.2        0.1280               3.5              3.00
## MSQ6_AG_399          4.6        0.1520               2.5              2.50
## MSQ6_AG_400          4.9        0.1290               3.0              4.00
## MSQ6_AG_401          5.1        0.1360               3.0              3.00
## MSQ6_AG_402          4.3        0.1480               1.5              2.00
## MSQ6_AG_403          4.9        0.1670               3.0              2.50
## MSQ6_AG_404          5.3        0.1180               2.5              3.00
## MSQ6_AG_405          5.8        0.1150               2.0              2.50
## MSQ6_AG_406          4.9        0.1240               3.0              3.00
## MSQ6_AG_407          5.2        0.1400               2.5              3.00
## MSQ6_AG_408          5.5        0.1320               3.0              4.00
## MSQ6_AG_409          5.7        0.1780               1.5              2.50
## MSQ6_AG_410          5.1        0.1140               1.0              2.00
## MSQ6_AG_411          5.4        0.1030               2.5              3.00
## MSQ6_AG_412          5.2        0.1430               4.0              4.00
## MSQ6_AG_413          4.9        0.0960               3.0              3.50
## MSQ6_AG_414          4.8        0.0820               1.5              1.00
## MSQ6_AG_415          4.5        0.0960               2.5              3.00
## MSQ6_AG_416          5.1        0.0860               2.0              2.50
## MSQ6_AG_417          4.8        0.1300               1.5              2.00
## MSQ6_AG_418          4.7        0.1200               3.0              4.00
## MSQ6_AG_419          4.2        0.1470               4.0              4.00
## MSQ6_AG_420          4.7        0.1230               3.5              3.50
## MSQ6_AG_421          4.3        0.1250               3.0              4.00
## MSQ6_AG_422          4.6        0.1420               2.5              3.00
## MSQ6_AG_423          4.4        0.1100               4.0              4.00
## MSQ6_AG_424          4.9        0.1830               3.5              4.00
## MSQ6_AG_425          4.3        0.1650               3.5              4.00
## MSQ6_AG_426          4.2        0.1450               4.0              3.00
## MSQ6_AG_427          7.8        0.0650               0.5              0.00
## MSQ6_AG_428          6.2        0.0930               0.0              0.00
## MSQ6_AG_429          6.9        0.0740               1.0              0.50
## MSQ6_AG_430          7.1        0.0710               0.0              0.00
## MSQ6_AG_431          6.8        0.0670               0.0              0.00
## MSQ6_AG_432          6.7        0.0820               0.0              0.00
## MSQ6_AG_433          6.5        0.0650               0.0              0.00
## MSQ6_AG_434          6.9        0.0650               0.0              0.50
## MSQ6_AG_490           NA            NA                NA                NA
## MSQ6_AG_491           NA            NA                NA                NA
## MSQ6_AG_492           NA            NA                NA                NA
## MSQ6_AG_493           NA            NA                NA                NA
## MSQ6_AG_494           NA            NA                NA                NA
## MSQ6_AG_495           NA            NA                NA                NA
## MSQ6_AG_496           NA            NA                NA                NA
## MSQ6_AG_497           NA            NA                NA                NA
## MSQ6_AG_498           NA            NA                NA                NA
## MSQ6_AG_499           NA            NA                NA                NA
## MSQ6_AG_500           NA            NA                NA                NA
## MSQ6_AG_501           NA            NA                NA                NA
## MSQ6_AG_502           NA            NA                NA                NA
## MSQ6_AG_503           NA            NA                NA                NA
## MSQ6_AG_504           NA            NA                NA                NA
## MSQ6_AG_505           NA            NA                NA                NA
## MSQ6_AG_506           NA            NA                NA                NA
## MSQ6_AG_509           NA            NA                NA                NA
## MSQ6_AG_511           NA            NA                NA                NA
## MSQ6_AG_512           NA            NA                NA                NA
## MSQ6_AG_514           NA            NA                NA                NA
## MSQ6_AG_515           NA            NA                NA                NA
## MSQ6_AG_516           NA            NA                NA                NA
## MSQ6_AG_517           NA            NA                NA                NA
## MSQ6_AG_518           NA            NA                NA                NA
## MSQ6_AG_519           NA            NA                NA                NA
## MSQ6_AG_520           NA            NA                NA                NA
## MSQ6_AG_521           NA            NA                NA                NA
## MSQ6_AG_522           NA            NA                NA                NA
## MSQ6_AG_523           NA            NA                NA                NA
## MSQ6_AG_524           NA            NA                NA                NA
## MSQ6_AG_525           NA            NA                NA                NA
## MSQ6_AG_526           NA            NA                NA                NA
## MSQ6_AG_527           NA            NA                NA                NA
## MSQ6_AG_528           NA            NA                NA                NA
## MSQ6_AG_529           NA            NA                NA                NA
## MSQ6_AG_530           NA            NA                NA                NA
##             histo_colo_total infiltr_CD3 infiltr_Gr1   CM_SCFA CM_BSCFA
## MSQ4_AG_012             1.00          NA          NA        NA       NA
## MSQ4_AG_013             0.00          NA          NA        NA       NA
## MSQ4_AG_014             0.00          NA          NA        NA       NA
## MSQ4_AG_015             0.00          NA          NA        NA       NA
## MSQ4_AG_016             1.00          NA          NA        NA       NA
## MSQ4_AG_017             0.00          NA          NA        NA       NA
## MSQ4_AG_018             0.00          NA          NA        NA       NA
## MSQ4_AG_019             1.00          NA          NA        NA       NA
## MSQ4_AG_020             0.00          NA          NA        NA       NA
## MSQ4_AG_021             0.00          NA          NA        NA       NA
## MSQ4_AG_022               NA          NA          NA        NA       NA
## MSQ4_AG_023             8.00          NA          NA        NA       NA
## MSQ4_AG_024             7.50          NA          NA        NA       NA
## MSQ4_AG_025             7.00          NA          NA        NA       NA
## MSQ4_AG_027             3.50          NA          NA        NA       NA
## MSQ4_AG_028             3.50          NA          NA        NA       NA
## MSQ4_AG_029             3.00          NA          NA        NA       NA
## MSQ4_AG_030             2.50          NA          NA        NA       NA
## MSQ4_AG_031             2.00          NA          NA        NA       NA
## MSQ4_AG_033             6.00          NA          NA        NA       NA
## MSQ4_AG_034             0.00          NA          NA        NA       NA
## MSQ4_AG_035             4.50          NA          NA        NA       NA
## MSQ4_AG_036             2.50          NA          NA        NA       NA
## MSQ4_AG_037             7.50          NA          NA        NA       NA
## MSQ4_AG_038             8.00          NA          NA        NA       NA
## MSQ4_AG_039             6.50          NA          NA        NA       NA
## MSQ4_AG_040             6.00          NA          NA        NA       NA
## MSQ4_AG_041             8.00          NA          NA        NA       NA
## MSQ4_AG_042             7.75          NA          NA        NA       NA
## MSQ4_AG_043             0.50          NA          NA        NA       NA
## MSQ4_AG_044             7.50          NA          NA        NA       NA
## MSQ4_AG_045             2.00          NA          NA        NA       NA
## MSQ4_AG_046             6.00          NA          NA        NA       NA
## MSQ4_AG_047             1.00          NA          NA        NA       NA
## MSQ4_AG_049             0.00          NA          NA        NA       NA
## MSQ4_AG_050             0.00          NA          NA        NA       NA
## MSQ4_AG_051             1.00          NA          NA        NA       NA
## MSQ4_AG_052             0.00          NA          NA        NA       NA
## MSQ4_AG_053             0.00          NA          NA        NA       NA
## MSQ4_AG_054             1.00          NA          NA        NA       NA
## MSQ4_AG_055             0.00          NA          NA        NA       NA
## MSQ4_AG_056             0.00          NA          NA        NA       NA
## MSQ4_AG_057               NA          NA          NA        NA       NA
## MSQ4_AG_058             8.00          NA          NA        NA       NA
## MSQ4_AG_059             7.50          NA          NA        NA       NA
## MSQ4_AG_062             3.50          NA          NA        NA       NA
## MSQ4_AG_064             3.00          NA          NA        NA       NA
## MSQ4_AG_065             2.50          NA          NA        NA       NA
## MSQ4_AG_072             7.50          NA          NA        NA       NA
## MSQ4_AG_076             8.00          NA          NA        NA       NA
## MSQ4_AG_077             7.75          NA          NA        NA       NA
## MSQ4_AG_080             2.00          NA          NA        NA       NA
## MSQ4_AG_083             0.00          NA          NA        NA       NA
## MSQ4_AG_084             0.00          NA          NA        NA       NA
## MSQ4_AG_085             0.00          NA          NA        NA       NA
## MSQ4_AG_086             1.00          NA          NA        NA       NA
## MSQ4_AG_087             0.00          NA          NA        NA       NA
## MSQ4_AG_088             0.00          NA          NA        NA       NA
## MSQ4_AG_089             1.00          NA          NA        NA       NA
## MSQ4_AG_091             0.00          NA          NA        NA       NA
## MSQ4_AG_092             8.00          NA          NA        NA       NA
## MSQ4_AG_093             7.50          NA          NA        NA       NA
## MSQ4_AG_095             6.50          NA          NA        NA       NA
## MSQ4_AG_096             3.50          NA          NA        NA       NA
## MSQ4_AG_101             4.00          NA          NA        NA       NA
## MSQ4_AG_103             0.00          NA          NA        NA       NA
## MSQ4_AG_104             4.50          NA          NA        NA       NA
## MSQ4_AG_105             2.50          NA          NA        NA       NA
## MSQ4_AG_106             7.50          NA          NA        NA       NA
## MSQ4_AG_107             8.00          NA          NA        NA       NA
## MSQ4_AG_108             6.50          NA          NA        NA       NA
## MSQ4_AG_109             6.00          NA          NA        NA       NA
## MSQ4_AG_111             7.75          NA          NA        NA       NA
## MSQ4_AG_112             0.50          NA          NA        NA       NA
## MSQ4_AG_113             7.50          NA          NA        NA       NA
## MSQ4_AG_117             0.00          NA          NA        NA       NA
## MSQ4_AG_119             0.00          NA          NA        NA       NA
## MSQ4_AG_121             0.00          NA          NA        NA       NA
## MSQ4_AG_125             0.00          NA          NA        NA       NA
## MSQ4_AG_127             7.50          NA          NA        NA       NA
## MSQ4_AG_128             7.00          NA          NA        NA       NA
## MSQ4_AG_129             6.50          NA          NA        NA       NA
## MSQ4_AG_130             3.50          NA          NA        NA       NA
## MSQ4_AG_131             3.50          NA          NA        NA       NA
## MSQ4_AG_132             3.00          NA          NA        NA       NA
## MSQ4_AG_133             2.50          NA          NA        NA       NA
## MSQ4_AG_135             4.00          NA          NA        NA       NA
## MSQ4_AG_136             6.00          NA          NA        NA       NA
## MSQ4_AG_141             8.00          NA          NA        NA       NA
## MSQ4_AG_143             6.00          NA          NA        NA       NA
## MSQ4_AG_144             8.00          NA          NA        NA       NA
## MSQ4_AG_145             7.75          NA          NA        NA       NA
## MSQ4_AG_148             2.00          NA          NA        NA       NA
## MSQ4_AG_149             6.00          NA          NA        NA       NA
## MSQ4_AG_151             0.00          NA          NA  94.84286        0
## MSQ4_AG_152             0.00          NA          NA  78.35625        0
## MSQ4_AG_153             0.00          NA          NA  59.55000        0
## MSQ4_AG_154             1.00          NA          NA   0.00000        0
## MSQ4_AG_155             0.00          NA          NA  82.88571        0
## MSQ4_AG_156             0.00          NA          NA  56.25625        0
## MSQ4_AG_157             1.00          NA          NA  58.89583        0
## MSQ4_AG_159             0.00          NA          NA  54.56429        0
## MSQ4_AG_160             8.00          NA          NA  67.29167        0
## MSQ4_AG_161             7.50          NA          NA  53.72500        0
## MSQ4_AG_162             7.00          NA          NA  64.27000        0
## MSQ4_AG_163             6.50          NA          NA  46.27353        0
## MSQ4_AG_164             3.50          NA          NA  96.28333        0
## MSQ4_AG_165             3.50          NA          NA  71.78750        0
## MSQ4_AG_167             2.50          NA          NA  66.14444        0
## MSQ4_AG_168             2.00          NA          NA  66.64231        0
## MSQ4_AG_169             4.00          NA          NA  76.37333        0
## MSQ4_AG_170             6.00          NA          NA  85.84167        0
## MSQ4_AG_171             0.00          NA          NA  81.42857        0
## MSQ4_AG_172             4.50          NA          NA  74.47083        0
## MSQ4_AG_173             2.50          NA          NA  79.26429        0
## MSQ4_AG_175             7.50          NA          NA 104.95625        0
## MSQ4_AG_176             8.00          NA          NA 123.90000        0
## MSQ4_AG_177             6.50          NA          NA  60.25000        0
## MSQ4_AG_178             6.00          NA          NA  64.34091        0
## MSQ4_AG_179             8.00          NA          NA  71.85000        0
## MSQ4_AG_180             0.50          NA          NA 103.80714        0
## MSQ4_AG_181             7.50          NA          NA  41.45000        0
## MSQ4_AG_183             6.00          NA          NA  70.75312        0
## MSQ4_AG_184               NA          NA          NA        NA       NA
## MSQ4_AG_185               NA          NA          NA        NA       NA
## MSQ4_AG_186               NA          NA          NA        NA       NA
## MSQ4_AG_187               NA          NA          NA        NA       NA
## MSQ4_AG_188               NA          NA          NA        NA       NA
## MSQ4_AG_189               NA          NA          NA        NA       NA
## MSQ4_AG_191               NA          NA          NA        NA       NA
## MSQ4_AG_192               NA          NA          NA        NA       NA
## MSQ4_AG_193               NA          NA          NA        NA       NA
## MSQ4_AG_194               NA          NA          NA        NA       NA
## MSQ4_AG_195               NA          NA          NA        NA       NA
## MSQ4_AG_196               NA          NA          NA        NA       NA
## MSQ4_AG_197               NA          NA          NA        NA       NA
## MSQ4_AG_198               NA          NA          NA        NA       NA
## MSQ4_AG_201               NA          NA          NA        NA       NA
## MSQ4_AG_205               NA          NA          NA        NA       NA
## MSQ4_AG_207               NA          NA          NA        NA       NA
## MSQ4_AG_208               NA          NA          NA        NA       NA
## MSQ4_AG_209               NA          NA          NA        NA       NA
## MSQ4_AG_210               NA          NA          NA        NA       NA
## MSQ4_AG_211               NA          NA          NA        NA       NA
## MSQ4_AG_212               NA          NA          NA        NA       NA
## MSQ4_AG_213               NA          NA          NA        NA       NA
## MSQ4_AG_215               NA          NA          NA        NA       NA
## MSQ4_AG_216               NA          NA          NA        NA       NA
## MSQ4_AG_217               NA          NA          NA        NA       NA
## MSQ4_AG_218               NA          NA          NA        NA       NA
## MSQ4_AG_219               NA          NA          NA        NA       NA
## MSQ4_AG_220               NA          NA          NA        NA       NA
## MSQ4_AG_221               NA          NA          NA        NA       NA
## MSQ4_AG_223               NA          NA          NA        NA       NA
## MSQ4_AG_224               NA          NA          NA        NA       NA
## MSQ4_AG_225               NA          NA          NA        NA       NA
## MSQ4_AG_226               NA          NA          NA        NA       NA
## MSQ4_AG_227               NA          NA          NA        NA       NA
## MSQ4_AG_228               NA          NA          NA        NA       NA
## MSQ4_AG_229               NA          NA          NA        NA       NA
## MSQ4_AG_230               NA          NA          NA        NA       NA
## MSQ4_AG_231               NA          NA          NA        NA       NA
## MSQ4_AG_232               NA          NA          NA        NA       NA
## MSQ4_AG_233               NA          NA          NA        NA       NA
## MSQ4_AG_234               NA          NA          NA        NA       NA
## MSQ4_AG_235               NA          NA          NA        NA       NA
## MSQ4_AG_236               NA          NA          NA        NA       NA
## MSQ4_AG_237               NA          NA          NA        NA       NA
## MSQ4_AG_238               NA          NA          NA        NA       NA
## MSQ4_AG_239               NA          NA          NA        NA       NA
## MSQ4_AG_240               NA          NA          NA        NA       NA
## MSQ4_AG_241               NA          NA          NA        NA       NA
## MSQ4_AG_242               NA          NA          NA        NA       NA
## MSQ4_AG_243               NA          NA          NA        NA       NA
## MSQ4_AG_244               NA          NA          NA        NA       NA
## MSQ4_AG_245               NA          NA          NA        NA       NA
## MSQ4_AG_247               NA          NA          NA        NA       NA
## MSQ4_AG_248               NA          NA          NA        NA       NA
## MSQ4_AG_249               NA          NA          NA        NA       NA
## MSQ6_AG_250               NA          NA          NA        NA       NA
## MSQ6_AG_251               NA          NA          NA        NA       NA
## MSQ6_AG_252               NA          NA          NA        NA       NA
## MSQ6_AG_253               NA          NA          NA        NA       NA
## MSQ6_AG_254               NA          NA          NA        NA       NA
## MSQ6_AG_255               NA          NA          NA        NA       NA
## MSQ6_AG_256               NA          NA          NA        NA       NA
## MSQ6_AG_257               NA          NA          NA        NA       NA
## MSQ6_AG_258               NA          NA          NA        NA       NA
## MSQ6_AG_259               NA          NA          NA        NA       NA
## MSQ6_AG_260               NA          NA          NA        NA       NA
## MSQ6_AG_261               NA          NA          NA        NA       NA
## MSQ6_AG_262               NA          NA          NA        NA       NA
## MSQ6_AG_263               NA          NA          NA        NA       NA
## MSQ6_AG_264               NA          NA          NA        NA       NA
## MSQ6_AG_265               NA          NA          NA        NA       NA
## MSQ6_AG_266               NA          NA          NA        NA       NA
## MSQ6_AG_267               NA          NA          NA        NA       NA
## MSQ6_AG_268               NA          NA          NA        NA       NA
## MSQ6_AG_269               NA          NA          NA        NA       NA
## MSQ6_AG_270               NA          NA          NA        NA       NA
## MSQ6_AG_271               NA          NA          NA        NA       NA
## MSQ6_AG_272               NA          NA          NA        NA       NA
## MSQ6_AG_273               NA          NA          NA        NA       NA
## MSQ6_AG_274               NA          NA          NA        NA       NA
## MSQ6_AG_275               NA          NA          NA        NA       NA
## MSQ6_AG_276               NA          NA          NA        NA       NA
## MSQ6_AG_277               NA          NA          NA        NA       NA
## MSQ6_AG_278               NA          NA          NA        NA       NA
## MSQ6_AG_279               NA          NA          NA        NA       NA
## MSQ6_AG_281               NA          NA          NA        NA       NA
## MSQ6_AG_282               NA          NA          NA        NA       NA
## MSQ6_AG_283               NA          NA          NA        NA       NA
## MSQ6_AG_284               NA          NA          NA        NA       NA
## MSQ6_AG_285               NA          NA          NA        NA       NA
## MSQ6_AG_286               NA          NA          NA        NA       NA
## MSQ6_AG_287               NA          NA          NA        NA       NA
## MSQ6_AG_288               NA          NA          NA        NA       NA
## MSQ6_AG_290               NA          NA          NA        NA       NA
## MSQ6_AG_291               NA          NA          NA        NA       NA
## MSQ6_AG_292               NA          NA          NA        NA       NA
## MSQ6_AG_293               NA          NA          NA        NA       NA
## MSQ6_AG_294               NA          NA          NA        NA       NA
## MSQ6_AG_295               NA          NA          NA        NA       NA
## MSQ6_AG_296               NA          NA          NA        NA       NA
## MSQ6_AG_298               NA          NA          NA        NA       NA
## MSQ6_AG_299               NA          NA          NA        NA       NA
## MSQ6_AG_300               NA          NA          NA        NA       NA
## MSQ6_AG_301               NA          NA          NA        NA       NA
## MSQ6_AG_302               NA          NA          NA        NA       NA
## MSQ6_AG_303               NA          NA          NA        NA       NA
## MSQ6_AG_304               NA          NA          NA        NA       NA
## MSQ6_AG_305               NA          NA          NA        NA       NA
## MSQ6_AG_306               NA          NA          NA        NA       NA
## MSQ6_AG_307               NA          NA          NA        NA       NA
## MSQ6_AG_308               NA          NA          NA        NA       NA
## MSQ6_AG_309               NA          NA          NA        NA       NA
## MSQ6_AG_310               NA          NA          NA        NA       NA
## MSQ6_AG_311               NA          NA          NA        NA       NA
## MSQ6_AG_312               NA          NA          NA        NA       NA
## MSQ6_AG_313               NA          NA          NA        NA       NA
## MSQ6_AG_314               NA          NA          NA        NA       NA
## MSQ6_AG_315               NA          NA          NA        NA       NA
## MSQ6_AG_316               NA          NA          NA        NA       NA
## MSQ6_AG_317               NA          NA          NA        NA       NA
## MSQ6_AG_318               NA          NA          NA        NA       NA
## MSQ6_AG_319               NA          NA          NA        NA       NA
## MSQ6_AG_320               NA          NA          NA        NA       NA
## MSQ6_AG_321               NA          NA          NA        NA       NA
## MSQ6_AG_322               NA          NA          NA        NA       NA
## MSQ6_AG_323               NA          NA          NA        NA       NA
## MSQ6_AG_324               NA          NA          NA        NA       NA
## MSQ6_AG_325               NA          NA          NA        NA       NA
## MSQ6_AG_326               NA          NA          NA        NA       NA
## MSQ6_AG_327               NA          NA          NA        NA       NA
## MSQ6_AG_328               NA          NA          NA        NA       NA
## MSQ6_AG_329               NA          NA          NA        NA       NA
## MSQ6_AG_330               NA          NA          NA        NA       NA
## MSQ6_AG_331               NA          NA          NA        NA       NA
## MSQ6_AG_332               NA          NA          NA        NA       NA
## MSQ6_AG_333               NA          NA          NA        NA       NA
## MSQ6_AG_334               NA          NA          NA        NA       NA
## MSQ6_AG_335               NA          NA          NA        NA       NA
## MSQ6_AG_336               NA          NA          NA        NA       NA
## MSQ6_AG_337               NA          NA          NA        NA       NA
## MSQ6_AG_339               NA          NA          NA        NA       NA
## MSQ6_AG_340               NA          NA          NA        NA       NA
## MSQ6_AG_341               NA          NA          NA        NA       NA
## MSQ6_AG_342               NA          NA          NA        NA       NA
## MSQ6_AG_343               NA          NA          NA        NA       NA
## MSQ6_AG_344             2.50          50          65        NA       NA
## MSQ6_AG_346             4.50          49          27        NA       NA
## MSQ6_AG_347             4.50          26          30        NA       NA
## MSQ6_AG_348             2.00          72          12        NA       NA
## MSQ6_AG_349             2.50          26          76        NA       NA
## MSQ6_AG_350             4.00          46          10        NA       NA
## MSQ6_AG_351             6.50         125         130        NA       NA
## MSQ6_AG_352             8.00          89          76        NA       NA
## MSQ6_AG_353             5.50          74          80        NA       NA
## MSQ6_AG_354             6.50         134         145        NA       NA
## MSQ6_AG_355             5.00         186         155        NA       NA
## MSQ6_AG_356             7.00          92         101        NA       NA
## MSQ6_AG_357             6.00         119          89        NA       NA
## MSQ6_AG_358             3.50          78         123        NA       NA
## MSQ6_AG_359             5.50          59          86        NA       NA
## MSQ6_AG_360             5.50          68          65        NA       NA
## MSQ6_AG_361             4.50          72          76        NA       NA
## MSQ6_AG_362             6.00          55          53        NA       NA
## MSQ6_AG_363             5.50          89          27        NA       NA
## MSQ6_AG_364             3.00          59          65        NA       NA
## MSQ6_AG_365             5.50          90          76        NA       NA
## MSQ6_AG_366             8.00          56          80        NA       NA
## MSQ6_AG_367             6.50          40          42        NA       NA
## MSQ6_AG_368             5.50          59          89        NA       NA
## MSQ6_AG_369             4.50          84         103        NA       NA
## MSQ6_AG_370             3.50          99          43        NA       NA
## MSQ6_AG_371             7.00          73          26        NA       NA
## MSQ6_AG_372             8.00         135          98        NA       NA
## MSQ6_AG_373             7.00         100         142        NA       NA
## MSQ6_AG_374             7.00         144          79        NA       NA
## MSQ6_AG_375             5.50         132          85        NA       NA
## MSQ6_AG_376             8.00         112          97        NA       NA
## MSQ6_AG_377             7.50          97          79        NA       NA
## MSQ6_AG_378             7.50         122          87        NA       NA
## MSQ6_AG_379             7.00         154         115        NA       NA
## MSQ6_AG_380             0.50          15           5        NA       NA
## MSQ6_AG_381             0.00          17           8        NA       NA
## MSQ6_AG_382             1.50          32          14        NA       NA
## MSQ6_AG_383             0.00          18           2        NA       NA
## MSQ6_AG_384             0.00          25          17        NA       NA
## MSQ6_AG_385             0.00          10          33        NA       NA
## MSQ6_AG_386             0.00           8          10        NA       NA
## MSQ6_AG_387             0.50          14           6        NA       NA
## MSQ6_AG_388             3.50          36          44  91.01250        0
## MSQ6_AG_389             2.50          50          65  60.40000        0
## MSQ6_AG_390             4.50          49          27        NA       NA
## MSQ6_AG_391             4.50          26          30  53.46250        0
## MSQ6_AG_392             2.00          72          12 112.22000        0
## MSQ6_AG_393             2.50          26          76  73.69167        0
## MSQ6_AG_394             4.00          46          10  71.39000        0
## MSQ6_AG_395             6.50         125         130        NA       NA
## MSQ6_AG_396             8.00          89          76        NA       NA
## MSQ6_AG_397             5.50          74          80        NA       NA
## MSQ6_AG_398             6.50         134         145        NA       NA
## MSQ6_AG_399             5.00         186         155        NA       NA
## MSQ6_AG_400             7.00          92         101        NA       NA
## MSQ6_AG_401             6.00         119          89  71.01250        0
## MSQ6_AG_402             3.50          78         123  51.11333        0
## MSQ6_AG_403             5.50          59          86        NA       NA
## MSQ6_AG_404             5.50          68          65        NA       NA
## MSQ6_AG_405             4.50          72          76        NA       NA
## MSQ6_AG_406             6.00          55          53  56.28333        0
## MSQ6_AG_407             5.50          89          27  53.45000        0
## MSQ6_AG_408             7.00          76          54        NA       NA
## MSQ6_AG_409             4.00          62          75        NA       NA
## MSQ6_AG_410             3.00          59          65        NA       NA
## MSQ6_AG_411             5.50          90          76        NA       NA
## MSQ6_AG_412             8.00          56          80  81.56875        0
## MSQ6_AG_413             6.50          40          42  47.17222        0
## MSQ6_AG_414             2.50          98          58        NA       NA
## MSQ6_AG_415             5.50          59          89  71.92000        0
## MSQ6_AG_416             4.50          84         103        NA       NA
## MSQ6_AG_417             3.50          99          43  47.13750        0
## MSQ6_AG_418             7.00          73          26  83.95000        0
## MSQ6_AG_419             8.00         135          98  46.05526        0
## MSQ6_AG_420             7.00         100         142  47.83750        0
## MSQ6_AG_421             7.00         144          79  58.51667        0
## MSQ6_AG_422             5.50         132          85  36.63333        0
## MSQ6_AG_423             8.00         112          97  45.75000        0
## MSQ6_AG_424             7.50          97          79        NA       NA
## MSQ6_AG_425             7.50         122          87        NA       NA
## MSQ6_AG_426             7.00         154         115        NA       NA
## MSQ6_AG_427             0.50          15           5        NA       NA
## MSQ6_AG_428             0.00          17           8        NA       NA
## MSQ6_AG_429             1.50          32          14        NA       NA
## MSQ6_AG_430             0.00          18           2        NA       NA
## MSQ6_AG_431             0.00          25          17        NA       NA
## MSQ6_AG_432             0.00          10          33  30.45000        0
## MSQ6_AG_433             0.00           8          10  43.15000        0
## MSQ6_AG_434             0.50          14           6        NA       NA
## MSQ6_AG_490               NA          NA          NA        NA       NA
## MSQ6_AG_491               NA          NA          NA        NA       NA
## MSQ6_AG_492               NA          NA          NA        NA       NA
## MSQ6_AG_493               NA          NA          NA        NA       NA
## MSQ6_AG_494               NA          NA          NA        NA       NA
## MSQ6_AG_495               NA          NA          NA        NA       NA
## MSQ6_AG_496               NA          NA          NA        NA       NA
## MSQ6_AG_497               NA          NA          NA        NA       NA
## MSQ6_AG_498               NA          NA          NA        NA       NA
## MSQ6_AG_499               NA          NA          NA        NA       NA
## MSQ6_AG_500               NA          NA          NA        NA       NA
## MSQ6_AG_501               NA          NA          NA        NA       NA
## MSQ6_AG_502               NA          NA          NA        NA       NA
## MSQ6_AG_503               NA          NA          NA        NA       NA
## MSQ6_AG_504               NA          NA          NA        NA       NA
## MSQ6_AG_505               NA          NA          NA        NA       NA
## MSQ6_AG_506               NA          NA          NA        NA       NA
## MSQ6_AG_509               NA          NA          NA        NA       NA
## MSQ6_AG_511               NA          NA          NA        NA       NA
## MSQ6_AG_512               NA          NA          NA        NA       NA
## MSQ6_AG_514               NA          NA          NA        NA       NA
## MSQ6_AG_515               NA          NA          NA        NA       NA
## MSQ6_AG_516               NA          NA          NA        NA       NA
## MSQ6_AG_517               NA          NA          NA        NA       NA
## MSQ6_AG_518               NA          NA          NA        NA       NA
## MSQ6_AG_519               NA          NA          NA        NA       NA
## MSQ6_AG_520               NA          NA          NA        NA       NA
## MSQ6_AG_521               NA          NA          NA        NA       NA
## MSQ6_AG_522               NA          NA          NA        NA       NA
## MSQ6_AG_523               NA          NA          NA        NA       NA
## MSQ6_AG_524               NA          NA          NA        NA       NA
## MSQ6_AG_525               NA          NA          NA        NA       NA
## MSQ6_AG_526               NA          NA          NA        NA       NA
## MSQ6_AG_527               NA          NA          NA        NA       NA
## MSQ6_AG_528               NA          NA          NA        NA       NA
## MSQ6_AG_529               NA          NA          NA        NA       NA
## MSQ6_AG_530               NA          NA          NA        NA       NA
##              CM_total CM_acetate CM_propioNAte CM_butyrate CM_succiNAte
## MSQ4_AG_012        NA         NA            NA          NA           NA
## MSQ4_AG_013        NA         NA            NA          NA           NA
## MSQ4_AG_014        NA         NA            NA          NA           NA
## MSQ4_AG_015        NA         NA            NA          NA           NA
## MSQ4_AG_016        NA         NA            NA          NA           NA
## MSQ4_AG_017        NA         NA            NA          NA           NA
## MSQ4_AG_018        NA         NA            NA          NA           NA
## MSQ4_AG_019        NA         NA            NA          NA           NA
## MSQ4_AG_020        NA         NA            NA          NA           NA
## MSQ4_AG_021        NA         NA            NA          NA           NA
## MSQ4_AG_022        NA         NA            NA          NA           NA
## MSQ4_AG_023        NA         NA            NA          NA           NA
## MSQ4_AG_024        NA         NA            NA          NA           NA
## MSQ4_AG_025        NA         NA            NA          NA           NA
## MSQ4_AG_027        NA         NA            NA          NA           NA
## MSQ4_AG_028        NA         NA            NA          NA           NA
## MSQ4_AG_029        NA         NA            NA          NA           NA
## MSQ4_AG_030        NA         NA            NA          NA           NA
## MSQ4_AG_031        NA         NA            NA          NA           NA
## MSQ4_AG_033        NA         NA            NA          NA           NA
## MSQ4_AG_034        NA         NA            NA          NA           NA
## MSQ4_AG_035        NA         NA            NA          NA           NA
## MSQ4_AG_036        NA         NA            NA          NA           NA
## MSQ4_AG_037        NA         NA            NA          NA           NA
## MSQ4_AG_038        NA         NA            NA          NA           NA
## MSQ4_AG_039        NA         NA            NA          NA           NA
## MSQ4_AG_040        NA         NA            NA          NA           NA
## MSQ4_AG_041        NA         NA            NA          NA           NA
## MSQ4_AG_042        NA         NA            NA          NA           NA
## MSQ4_AG_043        NA         NA            NA          NA           NA
## MSQ4_AG_044        NA         NA            NA          NA           NA
## MSQ4_AG_045        NA         NA            NA          NA           NA
## MSQ4_AG_046        NA         NA            NA          NA           NA
## MSQ4_AG_047        NA         NA            NA          NA           NA
## MSQ4_AG_049        NA         NA            NA          NA           NA
## MSQ4_AG_050        NA         NA            NA          NA           NA
## MSQ4_AG_051        NA         NA            NA          NA           NA
## MSQ4_AG_052        NA         NA            NA          NA           NA
## MSQ4_AG_053        NA         NA            NA          NA           NA
## MSQ4_AG_054        NA         NA            NA          NA           NA
## MSQ4_AG_055        NA         NA            NA          NA           NA
## MSQ4_AG_056        NA         NA            NA          NA           NA
## MSQ4_AG_057        NA         NA            NA          NA           NA
## MSQ4_AG_058        NA         NA            NA          NA           NA
## MSQ4_AG_059        NA         NA            NA          NA           NA
## MSQ4_AG_062        NA         NA            NA          NA           NA
## MSQ4_AG_064        NA         NA            NA          NA           NA
## MSQ4_AG_065        NA         NA            NA          NA           NA
## MSQ4_AG_072        NA         NA            NA          NA           NA
## MSQ4_AG_076        NA         NA            NA          NA           NA
## MSQ4_AG_077        NA         NA            NA          NA           NA
## MSQ4_AG_080        NA         NA            NA          NA           NA
## MSQ4_AG_083        NA         NA            NA          NA           NA
## MSQ4_AG_084        NA         NA            NA          NA           NA
## MSQ4_AG_085        NA         NA            NA          NA           NA
## MSQ4_AG_086        NA         NA            NA          NA           NA
## MSQ4_AG_087        NA         NA            NA          NA           NA
## MSQ4_AG_088        NA         NA            NA          NA           NA
## MSQ4_AG_089        NA         NA            NA          NA           NA
## MSQ4_AG_091        NA         NA            NA          NA           NA
## MSQ4_AG_092        NA         NA            NA          NA           NA
## MSQ4_AG_093        NA         NA            NA          NA           NA
## MSQ4_AG_095        NA         NA            NA          NA           NA
## MSQ4_AG_096        NA         NA            NA          NA           NA
## MSQ4_AG_101        NA         NA            NA          NA           NA
## MSQ4_AG_103        NA         NA            NA          NA           NA
## MSQ4_AG_104        NA         NA            NA          NA           NA
## MSQ4_AG_105        NA         NA            NA          NA           NA
## MSQ4_AG_106        NA         NA            NA          NA           NA
## MSQ4_AG_107        NA         NA            NA          NA           NA
## MSQ4_AG_108        NA         NA            NA          NA           NA
## MSQ4_AG_109        NA         NA            NA          NA           NA
## MSQ4_AG_111        NA         NA            NA          NA           NA
## MSQ4_AG_112        NA         NA            NA          NA           NA
## MSQ4_AG_113        NA         NA            NA          NA           NA
## MSQ4_AG_117        NA         NA            NA          NA           NA
## MSQ4_AG_119        NA         NA            NA          NA           NA
## MSQ4_AG_121        NA         NA            NA          NA           NA
## MSQ4_AG_125        NA         NA            NA          NA           NA
## MSQ4_AG_127        NA         NA            NA          NA           NA
## MSQ4_AG_128        NA         NA            NA          NA           NA
## MSQ4_AG_129        NA         NA            NA          NA           NA
## MSQ4_AG_130        NA         NA            NA          NA           NA
## MSQ4_AG_131        NA         NA            NA          NA           NA
## MSQ4_AG_132        NA         NA            NA          NA           NA
## MSQ4_AG_133        NA         NA            NA          NA           NA
## MSQ4_AG_135        NA         NA            NA          NA           NA
## MSQ4_AG_136        NA         NA            NA          NA           NA
## MSQ4_AG_141        NA         NA            NA          NA           NA
## MSQ4_AG_143        NA         NA            NA          NA           NA
## MSQ4_AG_144        NA         NA            NA          NA           NA
## MSQ4_AG_145        NA         NA            NA          NA           NA
## MSQ4_AG_148        NA         NA            NA          NA           NA
## MSQ4_AG_149        NA         NA            NA          NA           NA
## MSQ4_AG_151 102.86429   65.25714      7.764286   21.821429     4.585714
## MSQ4_AG_152  95.91250   57.77500     10.243750   10.337500    14.762500
## MSQ4_AG_153  71.95556   45.82222      7.727778    6.000000     7.327778
## MSQ4_AG_154        NA         NA            NA          NA           NA
## MSQ4_AG_155  99.63571   56.36429      8.600000   17.921429     4.692857
## MSQ4_AG_156  62.39688   34.06875      7.406250   14.781250     4.075000
## MSQ4_AG_157  63.52083   40.12083      6.158333   12.616667     2.200000
## MSQ4_AG_159  58.35714   36.75000      7.235714   10.578571     2.167857
## MSQ4_AG_160  80.92917   46.72083      9.712500   10.858333     8.041667
## MSQ4_AG_161  59.39722   35.98611      6.200000   11.538889     2.016667
## MSQ4_AG_162  66.27000   43.30000      9.495000   11.475000     2.000000
## MSQ4_AG_163  49.43235   32.61471      5.864706    7.794118     2.223529
## MSQ4_AG_164 100.92222   62.63333      8.455556   25.194444     2.222222
## MSQ4_AG_165  80.61667   49.08333     10.162500   12.541667     5.904167
## MSQ4_AG_167  76.87778   47.43889      9.316667    9.388889    17.971429
## MSQ4_AG_168  73.32308   49.60769      8.061538    8.973077     3.057692
## MSQ4_AG_169  90.04333   54.32333      7.730000   14.320000    13.670000
## MSQ4_AG_170 103.28333   60.06667     14.275000   11.500000    14.308333
## MSQ4_AG_171 100.94286   58.98571     10.107143   12.335714    14.514286
## MSQ4_AG_172  82.74583   52.88333      8.020833   13.566667     5.929167
## MSQ4_AG_173  85.44643   54.70714     10.346429   14.210714     3.889286
## MSQ4_AG_175 112.81250   68.30625     10.462500   26.187500     4.693750
## MSQ4_AG_176 133.90000   86.91250     13.912500   23.075000     4.387500
## MSQ4_AG_177  67.93929   42.65000      8.103571    9.496429     5.685714
## MSQ4_AG_178  75.01818   41.85455     10.672727   11.813636     6.863636
## MSQ4_AG_179  89.53500   51.28000     10.185000   10.385000    12.815000
## MSQ4_AG_180 107.05000   71.88571     13.564286   18.357143     3.242857
## MSQ4_AG_181  53.54000   26.71000     10.230000    4.510000     5.120000
## MSQ4_AG_183  78.49688   44.59063     11.134375   15.028125     5.212500
## MSQ4_AG_184        NA         NA            NA          NA           NA
## MSQ4_AG_185        NA         NA            NA          NA           NA
## MSQ4_AG_186        NA         NA            NA          NA           NA
## MSQ4_AG_187        NA         NA            NA          NA           NA
## MSQ4_AG_188        NA         NA            NA          NA           NA
## MSQ4_AG_189        NA         NA            NA          NA           NA
## MSQ4_AG_191        NA         NA            NA          NA           NA
## MSQ4_AG_192        NA         NA            NA          NA           NA
## MSQ4_AG_193        NA         NA            NA          NA           NA
## MSQ4_AG_194        NA         NA            NA          NA           NA
## MSQ4_AG_195        NA         NA            NA          NA           NA
## MSQ4_AG_196        NA         NA            NA          NA           NA
## MSQ4_AG_197        NA         NA            NA          NA           NA
## MSQ4_AG_198        NA         NA            NA          NA           NA
## MSQ4_AG_201        NA         NA            NA          NA           NA
## MSQ4_AG_205        NA         NA            NA          NA           NA
## MSQ4_AG_207        NA         NA            NA          NA           NA
## MSQ4_AG_208        NA         NA            NA          NA           NA
## MSQ4_AG_209        NA         NA            NA          NA           NA
## MSQ4_AG_210        NA         NA            NA          NA           NA
## MSQ4_AG_211        NA         NA            NA          NA           NA
## MSQ4_AG_212        NA         NA            NA          NA           NA
## MSQ4_AG_213        NA         NA            NA          NA           NA
## MSQ4_AG_215        NA         NA            NA          NA           NA
## MSQ4_AG_216        NA         NA            NA          NA           NA
## MSQ4_AG_217        NA         NA            NA          NA           NA
## MSQ4_AG_218        NA         NA            NA          NA           NA
## MSQ4_AG_219        NA         NA            NA          NA           NA
## MSQ4_AG_220        NA         NA            NA          NA           NA
## MSQ4_AG_221        NA         NA            NA          NA           NA
## MSQ4_AG_223        NA         NA            NA          NA           NA
## MSQ4_AG_224        NA         NA            NA          NA           NA
## MSQ4_AG_225        NA         NA            NA          NA           NA
## MSQ4_AG_226        NA         NA            NA          NA           NA
## MSQ4_AG_227        NA         NA            NA          NA           NA
## MSQ4_AG_228        NA         NA            NA          NA           NA
## MSQ4_AG_229        NA         NA            NA          NA           NA
## MSQ4_AG_230        NA         NA            NA          NA           NA
## MSQ4_AG_231        NA         NA            NA          NA           NA
## MSQ4_AG_232        NA         NA            NA          NA           NA
## MSQ4_AG_233        NA         NA            NA          NA           NA
## MSQ4_AG_234        NA         NA            NA          NA           NA
## MSQ4_AG_235        NA         NA            NA          NA           NA
## MSQ4_AG_236        NA         NA            NA          NA           NA
## MSQ4_AG_237        NA         NA            NA          NA           NA
## MSQ4_AG_238        NA         NA            NA          NA           NA
## MSQ4_AG_239        NA         NA            NA          NA           NA
## MSQ4_AG_240        NA         NA            NA          NA           NA
## MSQ4_AG_241        NA         NA            NA          NA           NA
## MSQ4_AG_242        NA         NA            NA          NA           NA
## MSQ4_AG_243        NA         NA            NA          NA           NA
## MSQ4_AG_244        NA         NA            NA          NA           NA
## MSQ4_AG_245        NA         NA            NA          NA           NA
## MSQ4_AG_247        NA         NA            NA          NA           NA
## MSQ4_AG_248        NA         NA            NA          NA           NA
## MSQ4_AG_249        NA         NA            NA          NA           NA
## MSQ6_AG_250        NA         NA            NA          NA           NA
## MSQ6_AG_251        NA         NA            NA          NA           NA
## MSQ6_AG_252        NA         NA            NA          NA           NA
## MSQ6_AG_253        NA         NA            NA          NA           NA
## MSQ6_AG_254        NA         NA            NA          NA           NA
## MSQ6_AG_255        NA         NA            NA          NA           NA
## MSQ6_AG_256        NA         NA            NA          NA           NA
## MSQ6_AG_257        NA         NA            NA          NA           NA
## MSQ6_AG_258        NA         NA            NA          NA           NA
## MSQ6_AG_259        NA         NA            NA          NA           NA
## MSQ6_AG_260        NA         NA            NA          NA           NA
## MSQ6_AG_261        NA         NA            NA          NA           NA
## MSQ6_AG_262        NA         NA            NA          NA           NA
## MSQ6_AG_263        NA         NA            NA          NA           NA
## MSQ6_AG_264        NA         NA            NA          NA           NA
## MSQ6_AG_265        NA         NA            NA          NA           NA
## MSQ6_AG_266        NA         NA            NA          NA           NA
## MSQ6_AG_267        NA         NA            NA          NA           NA
## MSQ6_AG_268        NA         NA            NA          NA           NA
## MSQ6_AG_269        NA         NA            NA          NA           NA
## MSQ6_AG_270        NA         NA            NA          NA           NA
## MSQ6_AG_271        NA         NA            NA          NA           NA
## MSQ6_AG_272        NA         NA            NA          NA           NA
## MSQ6_AG_273        NA         NA            NA          NA           NA
## MSQ6_AG_274        NA         NA            NA          NA           NA
## MSQ6_AG_275        NA         NA            NA          NA           NA
## MSQ6_AG_276        NA         NA            NA          NA           NA
## MSQ6_AG_277        NA         NA            NA          NA           NA
## MSQ6_AG_278        NA         NA            NA          NA           NA
## MSQ6_AG_279        NA         NA            NA          NA           NA
## MSQ6_AG_281        NA         NA            NA          NA           NA
## MSQ6_AG_282        NA         NA            NA          NA           NA
## MSQ6_AG_283        NA         NA            NA          NA           NA
## MSQ6_AG_284        NA         NA            NA          NA           NA
## MSQ6_AG_285        NA         NA            NA          NA           NA
## MSQ6_AG_286        NA         NA            NA          NA           NA
## MSQ6_AG_287        NA         NA            NA          NA           NA
## MSQ6_AG_288        NA         NA            NA          NA           NA
## MSQ6_AG_290        NA         NA            NA          NA           NA
## MSQ6_AG_291        NA         NA            NA          NA           NA
## MSQ6_AG_292        NA         NA            NA          NA           NA
## MSQ6_AG_293        NA         NA            NA          NA           NA
## MSQ6_AG_294        NA         NA            NA          NA           NA
## MSQ6_AG_295        NA         NA            NA          NA           NA
## MSQ6_AG_296        NA         NA            NA          NA           NA
## MSQ6_AG_298        NA         NA            NA          NA           NA
## MSQ6_AG_299        NA         NA            NA          NA           NA
## MSQ6_AG_300        NA         NA            NA          NA           NA
## MSQ6_AG_301        NA         NA            NA          NA           NA
## MSQ6_AG_302        NA         NA            NA          NA           NA
## MSQ6_AG_303        NA         NA            NA          NA           NA
## MSQ6_AG_304        NA         NA            NA          NA           NA
## MSQ6_AG_305        NA         NA            NA          NA           NA
## MSQ6_AG_306        NA         NA            NA          NA           NA
## MSQ6_AG_307        NA         NA            NA          NA           NA
## MSQ6_AG_308        NA         NA            NA          NA           NA
## MSQ6_AG_309        NA         NA            NA          NA           NA
## MSQ6_AG_310        NA         NA            NA          NA           NA
## MSQ6_AG_311        NA         NA            NA          NA           NA
## MSQ6_AG_312        NA         NA            NA          NA           NA
## MSQ6_AG_313        NA         NA            NA          NA           NA
## MSQ6_AG_314        NA         NA            NA          NA           NA
## MSQ6_AG_315        NA         NA            NA          NA           NA
## MSQ6_AG_316        NA         NA            NA          NA           NA
## MSQ6_AG_317        NA         NA            NA          NA           NA
## MSQ6_AG_318        NA         NA            NA          NA           NA
## MSQ6_AG_319        NA         NA            NA          NA           NA
## MSQ6_AG_320        NA         NA            NA          NA           NA
## MSQ6_AG_321        NA         NA            NA          NA           NA
## MSQ6_AG_322        NA         NA            NA          NA           NA
## MSQ6_AG_323        NA         NA            NA          NA           NA
## MSQ6_AG_324        NA         NA            NA          NA           NA
## MSQ6_AG_325        NA         NA            NA          NA           NA
## MSQ6_AG_326        NA         NA            NA          NA           NA
## MSQ6_AG_327        NA         NA            NA          NA           NA
## MSQ6_AG_328        NA         NA            NA          NA           NA
## MSQ6_AG_329        NA         NA            NA          NA           NA
## MSQ6_AG_330        NA         NA            NA          NA           NA
## MSQ6_AG_331        NA         NA            NA          NA           NA
## MSQ6_AG_332        NA         NA            NA          NA           NA
## MSQ6_AG_333        NA         NA            NA          NA           NA
## MSQ6_AG_334        NA         NA            NA          NA           NA
## MSQ6_AG_335        NA         NA            NA          NA           NA
## MSQ6_AG_336        NA         NA            NA          NA           NA
## MSQ6_AG_337        NA         NA            NA          NA           NA
## MSQ6_AG_339        NA         NA            NA          NA           NA
## MSQ6_AG_340        NA         NA            NA          NA           NA
## MSQ6_AG_341        NA         NA            NA          NA           NA
## MSQ6_AG_342        NA         NA            NA          NA           NA
## MSQ6_AG_343        NA         NA            NA          NA           NA
## MSQ6_AG_344        NA         NA            NA          NA           NA
## MSQ6_AG_346        NA         NA            NA          NA           NA
## MSQ6_AG_347        NA         NA            NA          NA           NA
## MSQ6_AG_348        NA         NA            NA          NA           NA
## MSQ6_AG_349        NA         NA            NA          NA           NA
## MSQ6_AG_350        NA         NA            NA          NA           NA
## MSQ6_AG_351        NA         NA            NA          NA           NA
## MSQ6_AG_352        NA         NA            NA          NA           NA
## MSQ6_AG_353        NA         NA            NA          NA           NA
## MSQ6_AG_354        NA         NA            NA          NA           NA
## MSQ6_AG_355        NA         NA            NA          NA           NA
## MSQ6_AG_356        NA         NA            NA          NA           NA
## MSQ6_AG_357        NA         NA            NA          NA           NA
## MSQ6_AG_358        NA         NA            NA          NA           NA
## MSQ6_AG_359        NA         NA            NA          NA           NA
## MSQ6_AG_360        NA         NA            NA          NA           NA
## MSQ6_AG_361        NA         NA            NA          NA           NA
## MSQ6_AG_362        NA         NA            NA          NA           NA
## MSQ6_AG_363        NA         NA            NA          NA           NA
## MSQ6_AG_364        NA         NA            NA          NA           NA
## MSQ6_AG_365        NA         NA            NA          NA           NA
## MSQ6_AG_366        NA         NA            NA          NA           NA
## MSQ6_AG_367        NA         NA            NA          NA           NA
## MSQ6_AG_368        NA         NA            NA          NA           NA
## MSQ6_AG_369        NA         NA            NA          NA           NA
## MSQ6_AG_370        NA         NA            NA          NA           NA
## MSQ6_AG_371        NA         NA            NA          NA           NA
## MSQ6_AG_372        NA         NA            NA          NA           NA
## MSQ6_AG_373        NA         NA            NA          NA           NA
## MSQ6_AG_374        NA         NA            NA          NA           NA
## MSQ6_AG_375        NA         NA            NA          NA           NA
## MSQ6_AG_376        NA         NA            NA          NA           NA
## MSQ6_AG_377        NA         NA            NA          NA           NA
## MSQ6_AG_378        NA         NA            NA          NA           NA
## MSQ6_AG_379        NA         NA            NA          NA           NA
## MSQ6_AG_380        NA         NA            NA          NA           NA
## MSQ6_AG_381        NA         NA            NA          NA           NA
## MSQ6_AG_382        NA         NA            NA          NA           NA
## MSQ6_AG_383        NA         NA            NA          NA           NA
## MSQ6_AG_384        NA         NA            NA          NA           NA
## MSQ6_AG_385        NA         NA            NA          NA           NA
## MSQ6_AG_386        NA         NA            NA          NA           NA
## MSQ6_AG_387        NA         NA            NA          NA           NA
## MSQ6_AG_388  97.96875   62.76875     16.281250   11.962500     3.018750
## MSQ6_AG_389  69.10833   41.15833     12.750000    6.491667     8.708333
## MSQ6_AG_390        NA         NA            NA          NA           NA
## MSQ6_AG_391  65.00000   34.08125     10.650000    8.731250     7.693750
## MSQ6_AG_392 114.99000   71.66000     19.720000   20.840000     0.000000
## MSQ6_AG_393  86.30000   49.57500     13.791667   10.325000     6.625000
## MSQ6_AG_394  78.31000   45.50000     16.790000    9.100000     3.500000
## MSQ6_AG_395        NA         NA            NA          NA           NA
## MSQ6_AG_396        NA         NA            NA          NA           NA
## MSQ6_AG_397        NA         NA            NA          NA           NA
## MSQ6_AG_398        NA         NA            NA          NA           NA
## MSQ6_AG_399        NA         NA            NA          NA           NA
## MSQ6_AG_400        NA         NA            NA          NA           NA
## MSQ6_AG_401  77.37500   47.56250      8.687500   14.762500     3.175000
## MSQ6_AG_402  60.86000   33.68333      7.146667   10.283333     7.033333
## MSQ6_AG_403        NA         NA            NA          NA           NA
## MSQ6_AG_404        NA         NA            NA          NA           NA
## MSQ6_AG_405        NA         NA            NA          NA           NA
## MSQ6_AG_406  56.28333   38.36667      9.383333    8.533333     0.000000
## MSQ6_AG_407  59.53333   34.16667     10.400000    8.883333     6.083333
## MSQ6_AG_408        NA         NA            NA          NA           NA
## MSQ6_AG_409        NA         NA            NA          NA           NA
## MSQ6_AG_410        NA         NA            NA          NA           NA
## MSQ6_AG_411        NA         NA            NA          NA           NA
## MSQ6_AG_412  86.06875   54.29375     14.206250   13.068750     2.156250
## MSQ6_AG_413  52.74444   31.54444      8.327778    7.300000     3.866667
## MSQ6_AG_414        NA         NA            NA          NA           NA
## MSQ6_AG_415  94.94000   55.75000      9.260000    6.910000     7.000000
## MSQ6_AG_416        NA         NA            NA          NA           NA
## MSQ6_AG_417  59.20833   32.78750      8.037500    6.312500     3.770833
## MSQ6_AG_418  90.26667   52.73333     14.491667   16.725000     3.650000
## MSQ6_AG_419  48.72368   31.61316      6.718421    6.568421     1.447368
## MSQ6_AG_420  50.44687   32.62812      5.525000    8.865625     1.350000
## MSQ6_AG_421  64.19167   36.34167     11.508333   10.666667     2.758333
## MSQ6_AG_422  45.25833   25.74167      6.316667    4.575000     3.908333
## MSQ6_AG_423  52.74091   26.85455      7.831818   11.063636     4.386364
## MSQ6_AG_424        NA         NA            NA          NA           NA
## MSQ6_AG_425        NA         NA            NA          NA           NA
## MSQ6_AG_426        NA         NA            NA          NA           NA
## MSQ6_AG_427        NA         NA            NA          NA           NA
## MSQ6_AG_428        NA         NA            NA          NA           NA
## MSQ6_AG_429        NA         NA            NA          NA           NA
## MSQ6_AG_430        NA         NA            NA          NA           NA
## MSQ6_AG_431        NA         NA            NA          NA           NA
## MSQ6_AG_432  35.41000   22.21000      4.370000    3.870000     4.960000
## MSQ6_AG_433  52.82500   33.08333      5.016667    5.050000     4.500000
## MSQ6_AG_434        NA         NA            NA          NA           NA
## MSQ6_AG_490        NA         NA            NA          NA           NA
## MSQ6_AG_491        NA         NA            NA          NA           NA
## MSQ6_AG_492        NA         NA            NA          NA           NA
## MSQ6_AG_493        NA         NA            NA          NA           NA
## MSQ6_AG_494        NA         NA            NA          NA           NA
## MSQ6_AG_495        NA         NA            NA          NA           NA
## MSQ6_AG_496        NA         NA            NA          NA           NA
## MSQ6_AG_497        NA         NA            NA          NA           NA
## MSQ6_AG_498        NA         NA            NA          NA           NA
## MSQ6_AG_499        NA         NA            NA          NA           NA
## MSQ6_AG_500        NA         NA            NA          NA           NA
## MSQ6_AG_501        NA         NA            NA          NA           NA
## MSQ6_AG_502        NA         NA            NA          NA           NA
## MSQ6_AG_503        NA         NA            NA          NA           NA
## MSQ6_AG_504        NA         NA            NA          NA           NA
## MSQ6_AG_505        NA         NA            NA          NA           NA
## MSQ6_AG_506        NA         NA            NA          NA           NA
## MSQ6_AG_509        NA         NA            NA          NA           NA
## MSQ6_AG_511        NA         NA            NA          NA           NA
## MSQ6_AG_512        NA         NA            NA          NA           NA
## MSQ6_AG_514        NA         NA            NA          NA           NA
## MSQ6_AG_515        NA         NA            NA          NA           NA
## MSQ6_AG_516        NA         NA            NA          NA           NA
## MSQ6_AG_517        NA         NA            NA          NA           NA
## MSQ6_AG_518        NA         NA            NA          NA           NA
## MSQ6_AG_519        NA         NA            NA          NA           NA
## MSQ6_AG_520        NA         NA            NA          NA           NA
## MSQ6_AG_521        NA         NA            NA          NA           NA
## MSQ6_AG_522        NA         NA            NA          NA           NA
## MSQ6_AG_523        NA         NA            NA          NA           NA
## MSQ6_AG_524        NA         NA            NA          NA           NA
## MSQ6_AG_525        NA         NA            NA          NA           NA
## MSQ6_AG_526        NA         NA            NA          NA           NA
## MSQ6_AG_527        NA         NA            NA          NA           NA
## MSQ6_AG_528        NA         NA            NA          NA           NA
## MSQ6_AG_529        NA         NA            NA          NA           NA
## MSQ6_AG_530        NA         NA            NA          NA           NA
##             CM_lactate   CB_total  CB_intact CB_permeable QPCR_total
## MSQ4_AG_012         NA         NA         NA           NA         NA
## MSQ4_AG_013         NA         NA         NA           NA         NA
## MSQ4_AG_014         NA         NA         NA           NA         NA
## MSQ4_AG_015         NA         NA         NA           NA         NA
## MSQ4_AG_016         NA         NA         NA           NA         NA
## MSQ4_AG_017         NA         NA         NA           NA         NA
## MSQ4_AG_018         NA         NA         NA           NA         NA
## MSQ4_AG_019         NA         NA         NA           NA         NA
## MSQ4_AG_020         NA         NA         NA           NA         NA
## MSQ4_AG_021         NA         NA         NA           NA         NA
## MSQ4_AG_022         NA         NA         NA           NA         NA
## MSQ4_AG_023         NA         NA         NA           NA         NA
## MSQ4_AG_024         NA         NA         NA           NA         NA
## MSQ4_AG_025         NA         NA         NA           NA         NA
## MSQ4_AG_027         NA         NA         NA           NA         NA
## MSQ4_AG_028         NA         NA         NA           NA         NA
## MSQ4_AG_029         NA         NA         NA           NA         NA
## MSQ4_AG_030         NA         NA         NA           NA         NA
## MSQ4_AG_031         NA         NA         NA           NA         NA
## MSQ4_AG_033         NA         NA         NA           NA         NA
## MSQ4_AG_034         NA         NA         NA           NA         NA
## MSQ4_AG_035         NA         NA         NA           NA         NA
## MSQ4_AG_036         NA         NA         NA           NA         NA
## MSQ4_AG_037         NA         NA         NA           NA         NA
## MSQ4_AG_038         NA         NA         NA           NA         NA
## MSQ4_AG_039         NA         NA         NA           NA         NA
## MSQ4_AG_040         NA         NA         NA           NA         NA
## MSQ4_AG_041         NA         NA         NA           NA         NA
## MSQ4_AG_042         NA         NA         NA           NA         NA
## MSQ4_AG_043         NA         NA         NA           NA         NA
## MSQ4_AG_044         NA         NA         NA           NA         NA
## MSQ4_AG_045         NA         NA         NA           NA         NA
## MSQ4_AG_046         NA         NA         NA           NA         NA
## MSQ4_AG_047         NA         NA         NA           NA         NA
## MSQ4_AG_049         NA         NA         NA           NA         NA
## MSQ4_AG_050         NA         NA         NA           NA         NA
## MSQ4_AG_051         NA         NA         NA           NA         NA
## MSQ4_AG_052         NA         NA         NA           NA         NA
## MSQ4_AG_053         NA         NA         NA           NA         NA
## MSQ4_AG_054         NA         NA         NA           NA         NA
## MSQ4_AG_055         NA         NA         NA           NA         NA
## MSQ4_AG_056         NA         NA         NA           NA         NA
## MSQ4_AG_057         NA         NA         NA           NA         NA
## MSQ4_AG_058         NA         NA         NA           NA         NA
## MSQ4_AG_059         NA         NA         NA           NA         NA
## MSQ4_AG_062         NA         NA         NA           NA         NA
## MSQ4_AG_064         NA         NA         NA           NA         NA
## MSQ4_AG_065         NA         NA         NA           NA         NA
## MSQ4_AG_072         NA         NA         NA           NA         NA
## MSQ4_AG_076         NA         NA         NA           NA         NA
## MSQ4_AG_077         NA         NA         NA           NA         NA
## MSQ4_AG_080         NA         NA         NA           NA         NA
## MSQ4_AG_083         NA         NA         NA           NA         NA
## MSQ4_AG_084         NA         NA         NA           NA         NA
## MSQ4_AG_085         NA         NA         NA           NA         NA
## MSQ4_AG_086         NA         NA         NA           NA         NA
## MSQ4_AG_087         NA         NA         NA           NA         NA
## MSQ4_AG_088         NA         NA         NA           NA         NA
## MSQ4_AG_089         NA         NA         NA           NA         NA
## MSQ4_AG_091         NA         NA         NA           NA         NA
## MSQ4_AG_092         NA         NA         NA           NA         NA
## MSQ4_AG_093         NA         NA         NA           NA         NA
## MSQ4_AG_095         NA         NA         NA           NA         NA
## MSQ4_AG_096         NA         NA         NA           NA         NA
## MSQ4_AG_101         NA         NA         NA           NA         NA
## MSQ4_AG_103         NA         NA         NA           NA         NA
## MSQ4_AG_104         NA         NA         NA           NA         NA
## MSQ4_AG_105         NA         NA         NA           NA         NA
## MSQ4_AG_106         NA         NA         NA           NA         NA
## MSQ4_AG_107         NA         NA         NA           NA         NA
## MSQ4_AG_108         NA         NA         NA           NA         NA
## MSQ4_AG_109         NA         NA         NA           NA         NA
## MSQ4_AG_111         NA         NA         NA           NA         NA
## MSQ4_AG_112         NA         NA         NA           NA         NA
## MSQ4_AG_113         NA         NA         NA           NA         NA
## MSQ4_AG_117         NA         NA         NA           NA         NA
## MSQ4_AG_119         NA         NA         NA           NA         NA
## MSQ4_AG_121         NA         NA         NA           NA         NA
## MSQ4_AG_125         NA         NA         NA           NA         NA
## MSQ4_AG_127         NA         NA         NA           NA         NA
## MSQ4_AG_128         NA         NA         NA           NA         NA
## MSQ4_AG_129         NA         NA         NA           NA         NA
## MSQ4_AG_130         NA         NA         NA           NA         NA
## MSQ4_AG_131         NA         NA         NA           NA         NA
## MSQ4_AG_132         NA         NA         NA           NA         NA
## MSQ4_AG_133         NA         NA         NA           NA         NA
## MSQ4_AG_135         NA         NA         NA           NA         NA
## MSQ4_AG_136         NA         NA         NA           NA         NA
## MSQ4_AG_141         NA         NA         NA           NA         NA
## MSQ4_AG_143         NA         NA         NA           NA         NA
## MSQ4_AG_144         NA         NA         NA           NA         NA
## MSQ4_AG_145         NA         NA         NA           NA         NA
## MSQ4_AG_148         NA         NA         NA           NA         NA
## MSQ4_AG_149         NA         NA         NA           NA         NA
## MSQ4_AG_151  3.4357143         NA         NA           NA   6.91e+08
## MSQ4_AG_152  2.7937500         NA         NA           NA   2.42e+08
## MSQ4_AG_153  5.0777778         NA         NA           NA   1.94e+09
## MSQ4_AG_154         NA         NA         NA           NA   4.76e+08
## MSQ4_AG_155 12.0571429         NA         NA           NA         NA
## MSQ4_AG_156  2.0656250         NA         NA           NA         NA
## MSQ4_AG_157  2.4250000         NA         NA           NA         NA
## MSQ4_AG_159  1.6250000         NA         NA           NA         NA
## MSQ4_AG_160  5.5958333         NA         NA           NA   9.17e+08
## MSQ4_AG_161  3.6555556         NA         NA           NA   1.23e+09
## MSQ4_AG_162  0.0000000         NA         NA           NA   1.19e+09
## MSQ4_AG_163  0.9352941         NA         NA           NA   9.21e+08
## MSQ4_AG_164  2.4166667         NA         NA           NA         NA
## MSQ4_AG_165  2.9250000         NA         NA           NA         NA
## MSQ4_AG_167  4.7035714         NA         NA           NA         NA
## MSQ4_AG_168  5.7333333         NA         NA           NA         NA
## MSQ4_AG_169  0.0000000         NA         NA           NA   5.75e+08
## MSQ4_AG_170  3.1333333         NA         NA           NA   7.12e+08
## MSQ4_AG_171  5.0000000         NA         NA           NA   6.89e+08
## MSQ4_AG_172  2.3458333         NA         NA           NA   6.42e+08
## MSQ4_AG_173  2.2928571         NA         NA           NA   4.17e+08
## MSQ4_AG_175  3.1625000         NA         NA           NA   8.95e+08
## MSQ4_AG_176  5.6125000         NA         NA           NA   5.65e+08
## MSQ4_AG_177  2.0035714         NA         NA           NA   7.92e+08
## MSQ4_AG_178  3.8136364         NA         NA           NA   9.34e+08
## MSQ4_AG_179  4.8700000         NA         NA           NA   7.90e+08
## MSQ4_AG_180  0.0000000         NA         NA           NA   1.16e+09
## MSQ4_AG_181  6.9700000         NA         NA           NA   6.63e+08
## MSQ4_AG_183  2.5312500         NA         NA           NA   8.27e+08
## MSQ4_AG_184         NA         NA         NA           NA         NA
## MSQ4_AG_185         NA         NA         NA           NA         NA
## MSQ4_AG_186         NA         NA         NA           NA         NA
## MSQ4_AG_187         NA         NA         NA           NA         NA
## MSQ4_AG_188         NA         NA         NA           NA         NA
## MSQ4_AG_189         NA         NA         NA           NA         NA
## MSQ4_AG_191         NA         NA         NA           NA         NA
## MSQ4_AG_192         NA         NA         NA           NA         NA
## MSQ4_AG_193         NA         NA         NA           NA         NA
## MSQ4_AG_194         NA         NA         NA           NA         NA
## MSQ4_AG_195         NA         NA         NA           NA         NA
## MSQ4_AG_196         NA         NA         NA           NA         NA
## MSQ4_AG_197         NA         NA         NA           NA         NA
## MSQ4_AG_198         NA         NA         NA           NA         NA
## MSQ4_AG_201         NA         NA         NA           NA         NA
## MSQ4_AG_205         NA         NA         NA           NA         NA
## MSQ4_AG_207         NA         NA         NA           NA         NA
## MSQ4_AG_208         NA         NA         NA           NA         NA
## MSQ4_AG_209         NA         NA         NA           NA         NA
## MSQ4_AG_210         NA         NA         NA           NA         NA
## MSQ4_AG_211         NA         NA         NA           NA         NA
## MSQ4_AG_212         NA         NA         NA           NA         NA
## MSQ4_AG_213         NA         NA         NA           NA         NA
## MSQ4_AG_215         NA         NA         NA           NA         NA
## MSQ4_AG_216         NA         NA         NA           NA         NA
## MSQ4_AG_217         NA         NA         NA           NA         NA
## MSQ4_AG_218         NA         NA         NA           NA         NA
## MSQ4_AG_219         NA         NA         NA           NA         NA
## MSQ4_AG_220         NA         NA         NA           NA         NA
## MSQ4_AG_221         NA         NA         NA           NA         NA
## MSQ4_AG_223         NA         NA         NA           NA         NA
## MSQ4_AG_224         NA         NA         NA           NA         NA
## MSQ4_AG_225         NA         NA         NA           NA         NA
## MSQ4_AG_226         NA         NA         NA           NA         NA
## MSQ4_AG_227         NA         NA         NA           NA         NA
## MSQ4_AG_228         NA         NA         NA           NA         NA
## MSQ4_AG_229         NA         NA         NA           NA         NA
## MSQ4_AG_230         NA         NA         NA           NA         NA
## MSQ4_AG_231         NA         NA         NA           NA         NA
## MSQ4_AG_232         NA         NA         NA           NA         NA
## MSQ4_AG_233         NA         NA         NA           NA         NA
## MSQ4_AG_234         NA         NA         NA           NA         NA
## MSQ4_AG_235         NA         NA         NA           NA         NA
## MSQ4_AG_236         NA         NA         NA           NA         NA
## MSQ4_AG_237         NA         NA         NA           NA         NA
## MSQ4_AG_238         NA         NA         NA           NA         NA
## MSQ4_AG_239         NA         NA         NA           NA         NA
## MSQ4_AG_240         NA         NA         NA           NA         NA
## MSQ4_AG_241         NA         NA         NA           NA         NA
## MSQ4_AG_242         NA         NA         NA           NA         NA
## MSQ4_AG_243         NA         NA         NA           NA         NA
## MSQ4_AG_244         NA         NA         NA           NA         NA
## MSQ4_AG_245         NA         NA         NA           NA         NA
## MSQ4_AG_247         NA         NA         NA           NA         NA
## MSQ4_AG_248         NA         NA         NA           NA         NA
## MSQ4_AG_249         NA         NA         NA           NA         NA
## MSQ6_AG_250         NA         NA         NA           NA         NA
## MSQ6_AG_251         NA         NA         NA           NA         NA
## MSQ6_AG_252         NA         NA         NA           NA         NA
## MSQ6_AG_253         NA         NA         NA           NA         NA
## MSQ6_AG_254         NA         NA         NA           NA         NA
## MSQ6_AG_255         NA         NA         NA           NA         NA
## MSQ6_AG_256         NA         NA         NA           NA         NA
## MSQ6_AG_257         NA         NA         NA           NA         NA
## MSQ6_AG_258         NA         NA         NA           NA         NA
## MSQ6_AG_259         NA         NA         NA           NA         NA
## MSQ6_AG_260         NA         NA         NA           NA         NA
## MSQ6_AG_261         NA         NA         NA           NA         NA
## MSQ6_AG_262         NA         NA         NA           NA         NA
## MSQ6_AG_263         NA         NA         NA           NA         NA
## MSQ6_AG_264         NA         NA         NA           NA         NA
## MSQ6_AG_265         NA         NA         NA           NA         NA
## MSQ6_AG_266         NA         NA         NA           NA         NA
## MSQ6_AG_267         NA         NA         NA           NA         NA
## MSQ6_AG_268         NA         NA         NA           NA         NA
## MSQ6_AG_269         NA         NA         NA           NA         NA
## MSQ6_AG_270         NA         NA         NA           NA         NA
## MSQ6_AG_271         NA         NA         NA           NA         NA
## MSQ6_AG_272         NA         NA         NA           NA         NA
## MSQ6_AG_273         NA         NA         NA           NA         NA
## MSQ6_AG_274         NA         NA         NA           NA         NA
## MSQ6_AG_275         NA         NA         NA           NA         NA
## MSQ6_AG_276         NA         NA         NA           NA         NA
## MSQ6_AG_277         NA         NA         NA           NA         NA
## MSQ6_AG_278         NA         NA         NA           NA         NA
## MSQ6_AG_279         NA         NA         NA           NA         NA
## MSQ6_AG_281         NA         NA         NA           NA         NA
## MSQ6_AG_282         NA         NA         NA           NA         NA
## MSQ6_AG_283         NA         NA         NA           NA         NA
## MSQ6_AG_284         NA         NA         NA           NA         NA
## MSQ6_AG_285         NA         NA         NA           NA         NA
## MSQ6_AG_286         NA         NA         NA           NA         NA
## MSQ6_AG_287         NA         NA         NA           NA         NA
## MSQ6_AG_288         NA         NA         NA           NA         NA
## MSQ6_AG_290         NA         NA         NA           NA         NA
## MSQ6_AG_291         NA         NA         NA           NA         NA
## MSQ6_AG_292         NA         NA         NA           NA         NA
## MSQ6_AG_293         NA         NA         NA           NA         NA
## MSQ6_AG_294         NA         NA         NA           NA         NA
## MSQ6_AG_295         NA         NA         NA           NA         NA
## MSQ6_AG_296         NA         NA         NA           NA         NA
## MSQ6_AG_298         NA         NA         NA           NA         NA
## MSQ6_AG_299         NA         NA         NA           NA         NA
## MSQ6_AG_300         NA         NA         NA           NA         NA
## MSQ6_AG_301         NA         NA         NA           NA         NA
## MSQ6_AG_302         NA         NA         NA           NA         NA
## MSQ6_AG_303         NA         NA         NA           NA         NA
## MSQ6_AG_304         NA         NA         NA           NA         NA
## MSQ6_AG_305         NA         NA         NA           NA         NA
## MSQ6_AG_306         NA         NA         NA           NA         NA
## MSQ6_AG_307         NA         NA         NA           NA         NA
## MSQ6_AG_308         NA         NA         NA           NA         NA
## MSQ6_AG_309         NA         NA         NA           NA         NA
## MSQ6_AG_310         NA         NA         NA           NA         NA
## MSQ6_AG_311         NA         NA         NA           NA         NA
## MSQ6_AG_312         NA         NA         NA           NA         NA
## MSQ6_AG_313         NA         NA         NA           NA         NA
## MSQ6_AG_314         NA         NA         NA           NA         NA
## MSQ6_AG_315         NA         NA         NA           NA         NA
## MSQ6_AG_316         NA         NA         NA           NA         NA
## MSQ6_AG_317         NA         NA         NA           NA         NA
## MSQ6_AG_318         NA         NA         NA           NA         NA
## MSQ6_AG_319         NA         NA         NA           NA         NA
## MSQ6_AG_320         NA         NA         NA           NA         NA
## MSQ6_AG_321         NA         NA         NA           NA         NA
## MSQ6_AG_322         NA         NA         NA           NA         NA
## MSQ6_AG_323         NA         NA         NA           NA         NA
## MSQ6_AG_324         NA         NA         NA           NA         NA
## MSQ6_AG_325         NA         NA         NA           NA         NA
## MSQ6_AG_326         NA         NA         NA           NA         NA
## MSQ6_AG_327         NA         NA         NA           NA         NA
## MSQ6_AG_328         NA         NA         NA           NA         NA
## MSQ6_AG_329         NA         NA         NA           NA         NA
## MSQ6_AG_330         NA         NA         NA           NA         NA
## MSQ6_AG_331         NA         NA         NA           NA         NA
## MSQ6_AG_332         NA         NA         NA           NA         NA
## MSQ6_AG_333         NA         NA         NA           NA         NA
## MSQ6_AG_334         NA         NA         NA           NA         NA
## MSQ6_AG_335         NA         NA         NA           NA         NA
## MSQ6_AG_336         NA         NA         NA           NA         NA
## MSQ6_AG_337         NA         NA         NA           NA         NA
## MSQ6_AG_339         NA         NA         NA           NA         NA
## MSQ6_AG_340         NA         NA         NA           NA         NA
## MSQ6_AG_341         NA         NA         NA           NA         NA
## MSQ6_AG_342         NA         NA         NA           NA         NA
## MSQ6_AG_343         NA         NA         NA           NA         NA
## MSQ6_AG_344         NA         NA         NA           NA         NA
## MSQ6_AG_346         NA         NA         NA           NA         NA
## MSQ6_AG_347         NA         NA         NA           NA         NA
## MSQ6_AG_348         NA         NA         NA           NA         NA
## MSQ6_AG_349         NA         NA         NA           NA         NA
## MSQ6_AG_350         NA         NA         NA           NA         NA
## MSQ6_AG_351         NA         NA         NA           NA         NA
## MSQ6_AG_352         NA         NA         NA           NA         NA
## MSQ6_AG_353         NA         NA         NA           NA         NA
## MSQ6_AG_354         NA         NA         NA           NA         NA
## MSQ6_AG_355         NA         NA         NA           NA         NA
## MSQ6_AG_356         NA         NA         NA           NA         NA
## MSQ6_AG_357         NA         NA         NA           NA         NA
## MSQ6_AG_358         NA         NA         NA           NA         NA
## MSQ6_AG_359         NA         NA         NA           NA         NA
## MSQ6_AG_360         NA         NA         NA           NA         NA
## MSQ6_AG_361         NA         NA         NA           NA         NA
## MSQ6_AG_362         NA         NA         NA           NA         NA
## MSQ6_AG_363         NA         NA         NA           NA         NA
## MSQ6_AG_364         NA         NA         NA           NA         NA
## MSQ6_AG_365         NA         NA         NA           NA         NA
## MSQ6_AG_366         NA         NA         NA           NA         NA
## MSQ6_AG_367         NA         NA         NA           NA         NA
## MSQ6_AG_368         NA         NA         NA           NA         NA
## MSQ6_AG_369         NA         NA         NA           NA         NA
## MSQ6_AG_370         NA         NA         NA           NA         NA
## MSQ6_AG_371         NA         NA         NA           NA         NA
## MSQ6_AG_372         NA         NA         NA           NA         NA
## MSQ6_AG_373         NA         NA         NA           NA         NA
## MSQ6_AG_374         NA         NA         NA           NA         NA
## MSQ6_AG_375         NA         NA         NA           NA         NA
## MSQ6_AG_376         NA         NA         NA           NA         NA
## MSQ6_AG_377         NA         NA         NA           NA         NA
## MSQ6_AG_378         NA         NA         NA           NA         NA
## MSQ6_AG_379         NA         NA         NA           NA         NA
## MSQ6_AG_380         NA         NA         NA           NA         NA
## MSQ6_AG_381         NA         NA         NA           NA         NA
## MSQ6_AG_382         NA         NA         NA           NA         NA
## MSQ6_AG_383         NA         NA         NA           NA         NA
## MSQ6_AG_384         NA         NA         NA           NA         NA
## MSQ6_AG_385         NA         NA         NA           NA         NA
## MSQ6_AG_386         NA         NA         NA           NA         NA
## MSQ6_AG_387         NA         NA         NA           NA         NA
## MSQ6_AG_388  3.9375000  561104439  310761183    234276594   7.38e+08
## MSQ6_AG_389  0.0000000  623292505  341053955    259174025   1.02e+09
## MSQ6_AG_390         NA  971157000  446151055    492746330   8.85e+08
## MSQ6_AG_391  3.8437500  480900780  239850407    226503780   9.83e+07
## MSQ6_AG_392  2.7700000 2165621920  805645280   1295483200   1.13e+09
## MSQ6_AG_393  5.9833333 1412034462  677288011    681150703   1.47e+09
## MSQ6_AG_394  3.4200000 1470938374  782951593    640588385   1.96e+09
## MSQ6_AG_395         NA 1139698637  733033956    358494231   1.22e+09
## MSQ6_AG_396         NA 1509446593  994474703    467501088   9.74e+08
## MSQ6_AG_397         NA 1676566833 1089115500    538974333   8.68e+08
## MSQ6_AG_398         NA 1269890165  760251462    481885099   1.62e+09
## MSQ6_AG_399         NA 1141410820  621695360    477942850   1.63e+09
## MSQ6_AG_400         NA 1076932758  544821176    479237604   1.28e+09
## MSQ6_AG_401  3.1875000 2065826714  896283506   1100427187   1.80e+09
## MSQ6_AG_402  2.7133333  851993389  466423467    354685100   8.95e+08
## MSQ6_AG_403         NA  454469400  217951631    223797623   2.44e+09
## MSQ6_AG_404         NA 1307821392  878586500    403043554   3.49e+09
## MSQ6_AG_405         NA  902581625  604359525    277697658   2.00e+09
## MSQ6_AG_406  0.0000000  751773275  482508583    248177142   1.25e+09
## MSQ6_AG_407  0.0000000  709642746  430873977    259873477   9.04e+08
## MSQ6_AG_408         NA 1105897414  534499114    538624114   1.92e+09
## MSQ6_AG_409         NA 1194567733  681376117    469278150   1.98e+09
## MSQ6_AG_410         NA 1242007433  682891978    504306978   2.04e+09
## MSQ6_AG_411         NA 1531467025  926221450    560804750   1.40e+08
## MSQ6_AG_412  2.3437500  811350430  490101040    293939360   7.35e+07
## MSQ6_AG_413  1.7055556 1115257963  688049588    375939575   9.40e+07
## MSQ6_AG_414         NA 1333631200  781537167    511545650   1.28e+08
## MSQ6_AG_415 16.0200000  888138075  443512850    394644800   9.87e+07
## MSQ6_AG_416         NA  195828967  103436667     81296967   2.03e+07
## MSQ6_AG_417  2.4291667  956652644  527511233    403194856   1.01e+08
## MSQ6_AG_418  2.6666667  940870975  496413042    394058683   9.85e+07
## MSQ6_AG_419  1.2210526  798202000  425939900    347543500   8.69e+07
## MSQ6_AG_420  1.2593750 1035975406  706501141    307601218   7.69e+07
## MSQ6_AG_421  2.9166667  935987580  533243480    371779980   9.29e+07
## MSQ6_AG_422  4.7166667 5617152750 3022237350   2325785550   5.81e+08
## MSQ6_AG_423  2.6045455  739727175  426222683    295072433   7.38e+07
## MSQ6_AG_424         NA 1018015500  596310200    394660300   9.87e+07
## MSQ6_AG_425         NA  819921300  339289683    349933467   8.75e+07
## MSQ6_AG_426         NA  123523950   98110650     22797775   5.70e+06
## MSQ6_AG_427         NA  766702444  428363589    310284822   7.76e+07
## MSQ6_AG_428         NA 1034705513  632954575    367522375   9.19e+07
## MSQ6_AG_429         NA 1147794120  682110330    426680210   1.07e+08
## MSQ6_AG_430         NA  820315980  454616580    337687130   8.44e+07
## MSQ6_AG_431         NA  964725300  391919733    531054333   1.33e+08
## MSQ6_AG_432  0.0000000 1164892300  571921167    549120183   1.37e+08
## MSQ6_AG_433  5.1750000 1331780817  800156500    482957017   1.21e+08
## MSQ6_AG_434         NA 1150284422  655875856    461378989   1.15e+08
## MSQ6_AG_490         NA         NA         NA           NA         NA
## MSQ6_AG_491         NA         NA         NA           NA         NA
## MSQ6_AG_492         NA         NA         NA           NA         NA
## MSQ6_AG_493         NA         NA         NA           NA         NA
## MSQ6_AG_494         NA         NA         NA           NA         NA
## MSQ6_AG_495         NA         NA         NA           NA         NA
## MSQ6_AG_496         NA         NA         NA           NA         NA
## MSQ6_AG_497         NA         NA         NA           NA         NA
## MSQ6_AG_498         NA         NA         NA           NA         NA
## MSQ6_AG_499         NA         NA         NA           NA         NA
## MSQ6_AG_500         NA         NA         NA           NA         NA
## MSQ6_AG_501         NA         NA         NA           NA         NA
## MSQ6_AG_502         NA         NA         NA           NA         NA
## MSQ6_AG_503         NA         NA         NA           NA         NA
## MSQ6_AG_504         NA         NA         NA           NA         NA
## MSQ6_AG_505         NA         NA         NA           NA         NA
## MSQ6_AG_506         NA         NA         NA           NA         NA
## MSQ6_AG_509         NA         NA         NA           NA         NA
## MSQ6_AG_511         NA         NA         NA           NA         NA
## MSQ6_AG_512         NA         NA         NA           NA         NA
## MSQ6_AG_514         NA         NA         NA           NA         NA
## MSQ6_AG_515         NA         NA         NA           NA         NA
## MSQ6_AG_516         NA         NA         NA           NA         NA
## MSQ6_AG_517         NA         NA         NA           NA         NA
## MSQ6_AG_518         NA         NA         NA           NA         NA
## MSQ6_AG_519         NA         NA         NA           NA         NA
## MSQ6_AG_520         NA         NA         NA           NA         NA
## MSQ6_AG_521         NA         NA         NA           NA         NA
## MSQ6_AG_522         NA         NA         NA           NA         NA
## MSQ6_AG_523         NA         NA         NA           NA         NA
## MSQ6_AG_524         NA         NA         NA           NA         NA
## MSQ6_AG_525         NA         NA         NA           NA         NA
## MSQ6_AG_526         NA         NA         NA           NA         NA
## MSQ6_AG_527         NA         NA         NA           NA         NA
## MSQ6_AG_528         NA         NA         NA           NA         NA
## MSQ6_AG_529         NA         NA         NA           NA         NA
## MSQ6_AG_530         NA         NA         NA           NA         NA
##             QPCR_Akkermansia QPCR_Enterobacteriaceae QPCR_Lachnospiraceae
## MSQ4_AG_012               NA                      NA                   NA
## MSQ4_AG_013               NA                      NA                   NA
## MSQ4_AG_014               NA                      NA                   NA
## MSQ4_AG_015               NA                      NA                   NA
## MSQ4_AG_016               NA                      NA                   NA
## MSQ4_AG_017               NA                      NA                   NA
## MSQ4_AG_018               NA                      NA                   NA
## MSQ4_AG_019               NA                      NA                   NA
## MSQ4_AG_020               NA                      NA                   NA
## MSQ4_AG_021               NA                      NA                   NA
## MSQ4_AG_022               NA                      NA                   NA
## MSQ4_AG_023               NA                      NA                   NA
## MSQ4_AG_024               NA                      NA                   NA
## MSQ4_AG_025               NA                      NA                   NA
## MSQ4_AG_027               NA                      NA                   NA
## MSQ4_AG_028               NA                      NA                   NA
## MSQ4_AG_029               NA                      NA                   NA
## MSQ4_AG_030               NA                      NA                   NA
## MSQ4_AG_031               NA                      NA                   NA
## MSQ4_AG_033               NA                      NA                   NA
## MSQ4_AG_034               NA                      NA                   NA
## MSQ4_AG_035               NA                      NA                   NA
## MSQ4_AG_036               NA                      NA                   NA
## MSQ4_AG_037               NA                      NA                   NA
## MSQ4_AG_038               NA                      NA                   NA
## MSQ4_AG_039               NA                      NA                   NA
## MSQ4_AG_040               NA                      NA                   NA
## MSQ4_AG_041               NA                      NA                   NA
## MSQ4_AG_042               NA                      NA                   NA
## MSQ4_AG_043               NA                      NA                   NA
## MSQ4_AG_044               NA                      NA                   NA
## MSQ4_AG_045               NA                      NA                   NA
## MSQ4_AG_046               NA                      NA                   NA
## MSQ4_AG_047               NA                      NA                   NA
## MSQ4_AG_049               NA                      NA                   NA
## MSQ4_AG_050               NA                      NA                   NA
## MSQ4_AG_051               NA                      NA                   NA
## MSQ4_AG_052               NA                      NA                   NA
## MSQ4_AG_053               NA                      NA                   NA
## MSQ4_AG_054               NA                      NA                   NA
## MSQ4_AG_055               NA                      NA                   NA
## MSQ4_AG_056               NA                      NA                   NA
## MSQ4_AG_057               NA                      NA                   NA
## MSQ4_AG_058               NA                      NA                   NA
## MSQ4_AG_059               NA                      NA                   NA
## MSQ4_AG_062               NA                      NA                   NA
## MSQ4_AG_064               NA                      NA                   NA
## MSQ4_AG_065               NA                      NA                   NA
## MSQ4_AG_072               NA                      NA                   NA
## MSQ4_AG_076               NA                      NA                   NA
## MSQ4_AG_077               NA                      NA                   NA
## MSQ4_AG_080               NA                      NA                   NA
## MSQ4_AG_083               NA                      NA                   NA
## MSQ4_AG_084               NA                      NA                   NA
## MSQ4_AG_085               NA                      NA                   NA
## MSQ4_AG_086               NA                      NA                   NA
## MSQ4_AG_087               NA                      NA                   NA
## MSQ4_AG_088               NA                      NA                   NA
## MSQ4_AG_089               NA                      NA                   NA
## MSQ4_AG_091               NA                      NA                   NA
## MSQ4_AG_092               NA                      NA                   NA
## MSQ4_AG_093               NA                      NA                   NA
## MSQ4_AG_095               NA                      NA                   NA
## MSQ4_AG_096               NA                      NA                   NA
## MSQ4_AG_101               NA                      NA                   NA
## MSQ4_AG_103               NA                      NA                   NA
## MSQ4_AG_104               NA                      NA                   NA
## MSQ4_AG_105               NA                      NA                   NA
## MSQ4_AG_106               NA                      NA                   NA
## MSQ4_AG_107               NA                      NA                   NA
## MSQ4_AG_108               NA                      NA                   NA
## MSQ4_AG_109               NA                      NA                   NA
## MSQ4_AG_111               NA                      NA                   NA
## MSQ4_AG_112               NA                      NA                   NA
## MSQ4_AG_113               NA                      NA                   NA
## MSQ4_AG_117               NA                      NA                   NA
## MSQ4_AG_119               NA                      NA                   NA
## MSQ4_AG_121               NA                      NA                   NA
## MSQ4_AG_125               NA                      NA                   NA
## MSQ4_AG_127               NA                      NA                   NA
## MSQ4_AG_128               NA                      NA                   NA
## MSQ4_AG_129               NA                      NA                   NA
## MSQ4_AG_130               NA                      NA                   NA
## MSQ4_AG_131               NA                      NA                   NA
## MSQ4_AG_132               NA                      NA                   NA
## MSQ4_AG_133               NA                      NA                   NA
## MSQ4_AG_135               NA                      NA                   NA
## MSQ4_AG_136               NA                      NA                   NA
## MSQ4_AG_141               NA                      NA                   NA
## MSQ4_AG_143               NA                      NA                   NA
## MSQ4_AG_144               NA                      NA                   NA
## MSQ4_AG_145               NA                      NA                   NA
## MSQ4_AG_148               NA                      NA                   NA
## MSQ4_AG_149               NA                      NA                   NA
## MSQ4_AG_151         2.42e+07                   20500             1.63e+08
## MSQ4_AG_152         3.17e+07                   15100             3.10e+07
## MSQ4_AG_153         5.84e+07                   32800             1.13e+08
## MSQ4_AG_154         1.80e+07                    4420             8.41e+07
## MSQ4_AG_155               NA                      NA                   NA
## MSQ4_AG_156               NA                      NA                   NA
## MSQ4_AG_157               NA                      NA                   NA
## MSQ4_AG_159               NA                      NA                   NA
## MSQ4_AG_160         2.25e+07                   27600             2.57e+07
## MSQ4_AG_161         2.87e+07                   11300             1.90e+08
## MSQ4_AG_162         5.32e+07                   31300             2.08e+08
## MSQ4_AG_163         3.17e+07                   24500             1.26e+08
## MSQ4_AG_164               NA                      NA                   NA
## MSQ4_AG_165               NA                      NA                   NA
## MSQ4_AG_167               NA                      NA                   NA
## MSQ4_AG_168               NA                      NA                   NA
## MSQ4_AG_169         1.70e+07                    9590             1.37e+08
## MSQ4_AG_170         2.36e+07                   12900             9.47e+07
## MSQ4_AG_171         6.76e+07                    6580             7.81e+07
## MSQ4_AG_172         1.74e+07                    3290             2.26e+08
## MSQ4_AG_173         1.25e+07                    3600             8.54e+07
## MSQ4_AG_175         3.36e+07                   27500             2.71e+08
## MSQ4_AG_176         1.81e+07                    7930             1.66e+08
## MSQ4_AG_177         3.61e+07                   13000             7.58e+07
## MSQ4_AG_178         5.22e+07                   25100             5.24e+07
## MSQ4_AG_179         5.73e+07                   26300             9.47e+07
## MSQ4_AG_180         4.21e+07                   80000             1.52e+08
## MSQ4_AG_181         3.84e+07                  167000             4.24e+07
## MSQ4_AG_183         2.69e+06                  257000             4.97e+07
## MSQ4_AG_184               NA                   32700                   NA
## MSQ4_AG_185               NA                      NA                   NA
## MSQ4_AG_186               NA                      NA                   NA
## MSQ4_AG_187               NA                      NA                   NA
## MSQ4_AG_188               NA                      NA                   NA
## MSQ4_AG_189               NA                      NA                   NA
## MSQ4_AG_191               NA                      NA                   NA
## MSQ4_AG_192               NA                      NA                   NA
## MSQ4_AG_193               NA                      NA                   NA
## MSQ4_AG_194               NA                      NA                   NA
## MSQ4_AG_195               NA                      NA                   NA
## MSQ4_AG_196               NA                      NA                   NA
## MSQ4_AG_197               NA                      NA                   NA
## MSQ4_AG_198               NA                      NA                   NA
## MSQ4_AG_201               NA                      NA                   NA
## MSQ4_AG_205               NA                      NA                   NA
## MSQ4_AG_207               NA                      NA                   NA
## MSQ4_AG_208               NA                      NA                   NA
## MSQ4_AG_209               NA                      NA                   NA
## MSQ4_AG_210               NA                      NA                   NA
## MSQ4_AG_211               NA                      NA                   NA
## MSQ4_AG_212               NA                      NA                   NA
## MSQ4_AG_213               NA                      NA                   NA
## MSQ4_AG_215               NA                      NA                   NA
## MSQ4_AG_216               NA                      NA                   NA
## MSQ4_AG_217               NA                      NA                   NA
## MSQ4_AG_218               NA                      NA                   NA
## MSQ4_AG_219               NA                      NA                   NA
## MSQ4_AG_220               NA                      NA                   NA
## MSQ4_AG_221               NA                      NA                   NA
## MSQ4_AG_223               NA                      NA                   NA
## MSQ4_AG_224               NA                      NA                   NA
## MSQ4_AG_225               NA                      NA                   NA
## MSQ4_AG_226               NA                      NA                   NA
## MSQ4_AG_227               NA                      NA                   NA
## MSQ4_AG_228               NA                      NA                   NA
## MSQ4_AG_229               NA                      NA                   NA
## MSQ4_AG_230               NA                      NA                   NA
## MSQ4_AG_231               NA                      NA                   NA
## MSQ4_AG_232               NA                      NA                   NA
## MSQ4_AG_233               NA                      NA                   NA
## MSQ4_AG_234               NA                      NA                   NA
## MSQ4_AG_235               NA                      NA                   NA
## MSQ4_AG_236               NA                      NA                   NA
## MSQ4_AG_237               NA                      NA                   NA
## MSQ4_AG_238               NA                      NA                   NA
## MSQ4_AG_239               NA                      NA                   NA
## MSQ4_AG_240               NA                      NA                   NA
## MSQ4_AG_241               NA                      NA                   NA
## MSQ4_AG_242               NA                      NA                   NA
## MSQ4_AG_243               NA                      NA                   NA
## MSQ4_AG_244               NA                      NA                   NA
## MSQ4_AG_245               NA                      NA                   NA
## MSQ4_AG_247               NA                      NA                   NA
## MSQ4_AG_248               NA                      NA                   NA
## MSQ4_AG_249               NA                      NA                   NA
## MSQ6_AG_250               NA                      NA                   NA
## MSQ6_AG_251               NA                      NA                   NA
## MSQ6_AG_252               NA                      NA                   NA
## MSQ6_AG_253               NA                      NA                   NA
## MSQ6_AG_254               NA                      NA                   NA
## MSQ6_AG_255               NA                      NA                   NA
## MSQ6_AG_256               NA                      NA                   NA
## MSQ6_AG_257               NA                      NA                   NA
## MSQ6_AG_258               NA                      NA                   NA
## MSQ6_AG_259               NA                      NA                   NA
## MSQ6_AG_260               NA                      NA                   NA
## MSQ6_AG_261               NA                      NA                   NA
## MSQ6_AG_262               NA                      NA                   NA
## MSQ6_AG_263               NA                      NA                   NA
## MSQ6_AG_264               NA                      NA                   NA
## MSQ6_AG_265               NA                      NA                   NA
## MSQ6_AG_266               NA                      NA                   NA
## MSQ6_AG_267               NA                      NA                   NA
## MSQ6_AG_268               NA                      NA                   NA
## MSQ6_AG_269               NA                      NA                   NA
## MSQ6_AG_270               NA                      NA                   NA
## MSQ6_AG_271               NA                      NA                   NA
## MSQ6_AG_272               NA                      NA                   NA
## MSQ6_AG_273               NA                      NA                   NA
## MSQ6_AG_274               NA                      NA                   NA
## MSQ6_AG_275               NA                      NA                   NA
## MSQ6_AG_276               NA                      NA                   NA
## MSQ6_AG_277               NA                      NA                   NA
## MSQ6_AG_278               NA                      NA                   NA
## MSQ6_AG_279               NA                      NA                   NA
## MSQ6_AG_281               NA                      NA                   NA
## MSQ6_AG_282               NA                      NA                   NA
## MSQ6_AG_283               NA                      NA                   NA
## MSQ6_AG_284               NA                      NA                   NA
## MSQ6_AG_285               NA                      NA                   NA
## MSQ6_AG_286               NA                      NA                   NA
## MSQ6_AG_287               NA                      NA                   NA
## MSQ6_AG_288               NA                      NA                   NA
## MSQ6_AG_290               NA                      NA                   NA
## MSQ6_AG_291               NA                      NA                   NA
## MSQ6_AG_292               NA                      NA                   NA
## MSQ6_AG_293               NA                      NA                   NA
## MSQ6_AG_294               NA                      NA                   NA
## MSQ6_AG_295               NA                      NA                   NA
## MSQ6_AG_296               NA                      NA                   NA
## MSQ6_AG_298               NA                      NA                   NA
## MSQ6_AG_299               NA                      NA                   NA
## MSQ6_AG_300               NA                      NA                   NA
## MSQ6_AG_301               NA                      NA                   NA
## MSQ6_AG_302               NA                      NA                   NA
## MSQ6_AG_303               NA                      NA                   NA
## MSQ6_AG_304               NA                      NA                   NA
## MSQ6_AG_305               NA                      NA                   NA
## MSQ6_AG_306               NA                      NA                   NA
## MSQ6_AG_307               NA                      NA                   NA
## MSQ6_AG_308               NA                      NA                   NA
## MSQ6_AG_309               NA                      NA                   NA
## MSQ6_AG_310               NA                      NA                   NA
## MSQ6_AG_311               NA                      NA                   NA
## MSQ6_AG_312               NA                      NA                   NA
## MSQ6_AG_313               NA                      NA                   NA
## MSQ6_AG_314               NA                      NA                   NA
## MSQ6_AG_315               NA                      NA                   NA
## MSQ6_AG_316               NA                      NA                   NA
## MSQ6_AG_317               NA                      NA                   NA
## MSQ6_AG_318               NA                      NA                   NA
## MSQ6_AG_319               NA                      NA                   NA
## MSQ6_AG_320               NA                      NA                   NA
## MSQ6_AG_321               NA                      NA                   NA
## MSQ6_AG_322               NA                      NA                   NA
## MSQ6_AG_323               NA                      NA                   NA
## MSQ6_AG_324               NA                      NA                   NA
## MSQ6_AG_325               NA                      NA                   NA
## MSQ6_AG_326               NA                      NA                   NA
## MSQ6_AG_327               NA                      NA                   NA
## MSQ6_AG_328               NA                      NA                   NA
## MSQ6_AG_329               NA                      NA                   NA
## MSQ6_AG_330               NA                      NA                   NA
## MSQ6_AG_331               NA                      NA                   NA
## MSQ6_AG_332               NA                      NA                   NA
## MSQ6_AG_333               NA                      NA                   NA
## MSQ6_AG_334               NA                      NA                   NA
## MSQ6_AG_335               NA                      NA                   NA
## MSQ6_AG_336               NA                      NA                   NA
## MSQ6_AG_337               NA                      NA                   NA
## MSQ6_AG_339               NA                      NA                   NA
## MSQ6_AG_340               NA                      NA                   NA
## MSQ6_AG_341               NA                      NA                   NA
## MSQ6_AG_342               NA                      NA                   NA
## MSQ6_AG_343               NA                      NA                   NA
## MSQ6_AG_344               NA                      NA                   NA
## MSQ6_AG_346               NA                      NA                   NA
## MSQ6_AG_347               NA                      NA                   NA
## MSQ6_AG_348               NA                      NA                   NA
## MSQ6_AG_349               NA                      NA                   NA
## MSQ6_AG_350               NA                      NA                   NA
## MSQ6_AG_351               NA                      NA                   NA
## MSQ6_AG_352               NA                      NA                   NA
## MSQ6_AG_353               NA                      NA                   NA
## MSQ6_AG_354               NA                      NA                   NA
## MSQ6_AG_355               NA                      NA                   NA
## MSQ6_AG_356               NA                      NA                   NA
## MSQ6_AG_357               NA                      NA                   NA
## MSQ6_AG_358               NA                      NA                   NA
## MSQ6_AG_359               NA                      NA                   NA
## MSQ6_AG_360               NA                      NA                   NA
## MSQ6_AG_361               NA                      NA                   NA
## MSQ6_AG_362               NA                      NA                   NA
## MSQ6_AG_363               NA                      NA                   NA
## MSQ6_AG_364               NA                      NA                   NA
## MSQ6_AG_365               NA                      NA                   NA
## MSQ6_AG_366               NA                      NA                   NA
## MSQ6_AG_367               NA                      NA                   NA
## MSQ6_AG_368               NA                      NA                   NA
## MSQ6_AG_369               NA                      NA                   NA
## MSQ6_AG_370               NA                      NA                   NA
## MSQ6_AG_371               NA                      NA                   NA
## MSQ6_AG_372               NA                      NA                   NA
## MSQ6_AG_373               NA                      NA                   NA
## MSQ6_AG_374               NA                      NA                   NA
## MSQ6_AG_375               NA                      NA                   NA
## MSQ6_AG_376               NA                      NA                   NA
## MSQ6_AG_377               NA                      NA                   NA
## MSQ6_AG_378               NA                      NA                   NA
## MSQ6_AG_379               NA                      NA                   NA
## MSQ6_AG_380               NA                      NA                   NA
## MSQ6_AG_381               NA                      NA                   NA
## MSQ6_AG_382               NA                      NA                   NA
## MSQ6_AG_383               NA                      NA                   NA
## MSQ6_AG_384               NA                      NA                   NA
## MSQ6_AG_385               NA                      NA                   NA
## MSQ6_AG_386               NA                      NA                   NA
## MSQ6_AG_387               NA                      NA                   NA
## MSQ6_AG_388         1.94e+07                   46900             1.01e+08
## MSQ6_AG_389         3.10e+07                  129000             8.31e+07
## MSQ6_AG_390         2.78e+07                    7950             1.43e+08
## MSQ6_AG_391         5.05e+06                    3350             3.25e+06
## MSQ6_AG_392         3.38e+07                    9790             1.66e+08
## MSQ6_AG_393         4.64e+07                   62200             1.08e+08
## MSQ6_AG_394         6.47e+07                   79000             2.85e+08
## MSQ6_AG_395         6.46e+07                   71200             1.15e+08
## MSQ6_AG_396         3.60e+07                   50800             9.31e+07
## MSQ6_AG_397         2.25e+07                  176000             5.57e+07
## MSQ6_AG_398         2.15e+07                   69400             2.90e+08
## MSQ6_AG_399         1.14e+08                  237000             3.04e+08
## MSQ6_AG_400         4.77e+07                  115000             9.57e+07
## MSQ6_AG_401         3.40e+07                   61300             3.20e+08
## MSQ6_AG_402         5.10e+07                   26000             6.65e+07
## MSQ6_AG_403         5.53e+07                  130000             2.45e+08
## MSQ6_AG_404         9.38e+07                 1190000             4.39e+08
## MSQ6_AG_405         7.72e+07                  539000             1.98e+08
## MSQ6_AG_406         6.60e+07                  973000             9.77e+07
## MSQ6_AG_407         3.37e+07                  231000             6.27e+07
## MSQ6_AG_408         5.63e+07                   80600             1.02e+08
## MSQ6_AG_409         6.97e+07                  156000             1.77e+08
## MSQ6_AG_410         6.88e+07                  537000             1.59e+08
## MSQ6_AG_411         5.71e+06                 1320000             2.44e+08
## MSQ6_AG_412         5.14e+07                 2770000             9.60e+07
## MSQ6_AG_413         7.56e+07                 1090000             1.41e+08
## MSQ6_AG_414         3.25e+07                   49600             1.72e+08
## MSQ6_AG_415         1.09e+07                  689000             1.48e+08
## MSQ6_AG_416         3.43e+06                 5210000             4.46e+07
## MSQ6_AG_417         3.08e+07                  358000             1.33e+08
## MSQ6_AG_418         4.52e+07                  333000             1.19e+08
## MSQ6_AG_419         2.38e+07                   85200             1.57e+08
## MSQ6_AG_420         4.06e+07                   23300             2.74e+08
## MSQ6_AG_421         5.92e+07                  383000             1.72e+08
## MSQ6_AG_422         1.53e+07                   12500             8.27e+07
## MSQ6_AG_423         7.26e+07                   65100             2.48e+08
## MSQ6_AG_424         3.34e+07                  352000             2.26e+08
## MSQ6_AG_425         4.61e+07                  199000             2.28e+08
## MSQ6_AG_426         1.05e+08                88600000             4.16e+08
## MSQ6_AG_427         3.17e+07                   63600             2.81e+08
## MSQ6_AG_428         3.56e+06                   36800             2.51e+08
## MSQ6_AG_429         1.31e+07                   56100             6.81e+08
## MSQ6_AG_430         1.31e+07                   54700             2.79e+08
## MSQ6_AG_431         1.83e+07                  280000             1.36e+09
## MSQ6_AG_432         6.58e+06                   48500             2.85e+08
## MSQ6_AG_433         6.22e+06                   16300             6.07e+07
## MSQ6_AG_434         3.68e+07                   73200             3.43e+08
## MSQ6_AG_490               NA                      NA                   NA
## MSQ6_AG_491               NA                      NA                   NA
## MSQ6_AG_492               NA                      NA                   NA
## MSQ6_AG_493               NA                      NA                   NA
## MSQ6_AG_494               NA                      NA                   NA
## MSQ6_AG_495               NA                      NA                   NA
## MSQ6_AG_496               NA                      NA                   NA
## MSQ6_AG_497               NA                      NA                   NA
## MSQ6_AG_498               NA                      NA                   NA
## MSQ6_AG_499               NA                      NA                   NA
## MSQ6_AG_500               NA                      NA                   NA
## MSQ6_AG_501               NA                      NA                   NA
## MSQ6_AG_502               NA                      NA                   NA
## MSQ6_AG_503               NA                      NA                   NA
## MSQ6_AG_504               NA                      NA                   NA
## MSQ6_AG_505               NA                      NA                   NA
## MSQ6_AG_506               NA                      NA                   NA
## MSQ6_AG_509               NA                      NA                   NA
## MSQ6_AG_511               NA                      NA                   NA
## MSQ6_AG_512               NA                      NA                   NA
## MSQ6_AG_514               NA                      NA                   NA
## MSQ6_AG_515               NA                      NA                   NA
## MSQ6_AG_516               NA                      NA                   NA
## MSQ6_AG_517               NA                      NA                   NA
## MSQ6_AG_518               NA                      NA                   NA
## MSQ6_AG_519               NA                      NA                   NA
## MSQ6_AG_520               NA                      NA                   NA
## MSQ6_AG_521               NA                      NA                   NA
## MSQ6_AG_522               NA                      NA                   NA
## MSQ6_AG_523               NA                      NA                   NA
## MSQ6_AG_524               NA                      NA                   NA
## MSQ6_AG_525               NA                      NA                   NA
## MSQ6_AG_526               NA                      NA                   NA
## MSQ6_AG_527               NA                      NA                   NA
## MSQ6_AG_528               NA                      NA                   NA
## MSQ6_AG_529               NA                      NA                   NA
## MSQ6_AG_530               NA                      NA                   NA
##             QPCR_Ruminococcaceae QPCR_ButCoA Tyramine Putrescine Cadaverine
## MSQ4_AG_012                   NA          NA       NA         NA         NA
## MSQ4_AG_013                   NA          NA       NA         NA         NA
## MSQ4_AG_014                   NA          NA       NA         NA         NA
## MSQ4_AG_015                   NA          NA       NA         NA         NA
## MSQ4_AG_016                   NA          NA       NA         NA         NA
## MSQ4_AG_017                   NA          NA       NA         NA         NA
## MSQ4_AG_018                   NA          NA       NA         NA         NA
## MSQ4_AG_019                   NA          NA       NA         NA         NA
## MSQ4_AG_020                   NA          NA       NA         NA         NA
## MSQ4_AG_021                   NA          NA       NA         NA         NA
## MSQ4_AG_022                   NA          NA       NA         NA         NA
## MSQ4_AG_023                   NA          NA       NA         NA         NA
## MSQ4_AG_024                   NA          NA       NA         NA         NA
## MSQ4_AG_025                   NA          NA       NA         NA         NA
## MSQ4_AG_027                   NA          NA       NA         NA         NA
## MSQ4_AG_028                   NA          NA       NA         NA         NA
## MSQ4_AG_029                   NA          NA       NA         NA         NA
## MSQ4_AG_030                   NA          NA       NA         NA         NA
## MSQ4_AG_031                   NA          NA       NA         NA         NA
## MSQ4_AG_033                   NA          NA       NA         NA         NA
## MSQ4_AG_034                   NA          NA       NA         NA         NA
## MSQ4_AG_035                   NA          NA       NA         NA         NA
## MSQ4_AG_036                   NA          NA       NA         NA         NA
## MSQ4_AG_037                   NA          NA       NA         NA         NA
## MSQ4_AG_038                   NA          NA       NA         NA         NA
## MSQ4_AG_039                   NA          NA       NA         NA         NA
## MSQ4_AG_040                   NA          NA       NA         NA         NA
## MSQ4_AG_041                   NA          NA       NA         NA         NA
## MSQ4_AG_042                   NA          NA       NA         NA         NA
## MSQ4_AG_043                   NA          NA       NA         NA         NA
## MSQ4_AG_044                   NA          NA       NA         NA         NA
## MSQ4_AG_045                   NA          NA       NA         NA         NA
## MSQ4_AG_046                   NA          NA       NA         NA         NA
## MSQ4_AG_047                   NA          NA       NA         NA         NA
## MSQ4_AG_049                   NA          NA       NA         NA         NA
## MSQ4_AG_050                   NA          NA       NA         NA         NA
## MSQ4_AG_051                   NA          NA       NA         NA         NA
## MSQ4_AG_052                   NA          NA       NA         NA         NA
## MSQ4_AG_053                   NA          NA       NA         NA         NA
## MSQ4_AG_054                   NA          NA       NA         NA         NA
## MSQ4_AG_055                   NA          NA       NA         NA         NA
## MSQ4_AG_056                   NA          NA       NA         NA         NA
## MSQ4_AG_057                   NA          NA       NA         NA         NA
## MSQ4_AG_058                   NA          NA       NA         NA         NA
## MSQ4_AG_059                   NA          NA       NA         NA         NA
## MSQ4_AG_062                   NA          NA       NA         NA         NA
## MSQ4_AG_064                   NA          NA       NA         NA         NA
## MSQ4_AG_065                   NA          NA       NA         NA         NA
## MSQ4_AG_072                   NA          NA       NA         NA         NA
## MSQ4_AG_076                   NA          NA       NA         NA         NA
## MSQ4_AG_077                   NA          NA       NA         NA         NA
## MSQ4_AG_080                   NA          NA       NA         NA         NA
## MSQ4_AG_083                   NA          NA       NA         NA         NA
## MSQ4_AG_084                   NA          NA       NA         NA         NA
## MSQ4_AG_085                   NA          NA       NA         NA         NA
## MSQ4_AG_086                   NA          NA       NA         NA         NA
## MSQ4_AG_087                   NA          NA       NA         NA         NA
## MSQ4_AG_088                   NA          NA       NA         NA         NA
## MSQ4_AG_089                   NA          NA       NA         NA         NA
## MSQ4_AG_091                   NA          NA       NA         NA         NA
## MSQ4_AG_092                   NA          NA       NA         NA         NA
## MSQ4_AG_093                   NA          NA       NA         NA         NA
## MSQ4_AG_095                   NA          NA       NA         NA         NA
## MSQ4_AG_096                   NA          NA       NA         NA         NA
## MSQ4_AG_101                   NA          NA       NA         NA         NA
## MSQ4_AG_103                   NA          NA       NA         NA         NA
## MSQ4_AG_104                   NA          NA       NA         NA         NA
## MSQ4_AG_105                   NA          NA       NA         NA         NA
## MSQ4_AG_106                   NA          NA       NA         NA         NA
## MSQ4_AG_107                   NA          NA       NA         NA         NA
## MSQ4_AG_108                   NA          NA       NA         NA         NA
## MSQ4_AG_109                   NA          NA       NA         NA         NA
## MSQ4_AG_111                   NA          NA       NA         NA         NA
## MSQ4_AG_112                   NA          NA       NA         NA         NA
## MSQ4_AG_113                   NA          NA       NA         NA         NA
## MSQ4_AG_117                   NA          NA       NA         NA         NA
## MSQ4_AG_119                   NA          NA       NA         NA         NA
## MSQ4_AG_121                   NA          NA       NA         NA         NA
## MSQ4_AG_125                   NA          NA       NA         NA         NA
## MSQ4_AG_127                   NA          NA       NA         NA         NA
## MSQ4_AG_128                   NA          NA       NA         NA         NA
## MSQ4_AG_129                   NA          NA       NA         NA         NA
## MSQ4_AG_130                   NA          NA       NA         NA         NA
## MSQ4_AG_131                   NA          NA       NA         NA         NA
## MSQ4_AG_132                   NA          NA       NA         NA         NA
## MSQ4_AG_133                   NA          NA       NA         NA         NA
## MSQ4_AG_135                   NA          NA       NA         NA         NA
## MSQ4_AG_136                   NA          NA       NA         NA         NA
## MSQ4_AG_141                   NA          NA       NA         NA         NA
## MSQ4_AG_143                   NA          NA       NA         NA         NA
## MSQ4_AG_144                   NA          NA       NA         NA         NA
## MSQ4_AG_145                   NA          NA       NA         NA         NA
## MSQ4_AG_148                   NA          NA       NA         NA         NA
## MSQ4_AG_149                   NA          NA       NA         NA         NA
## MSQ4_AG_151             2.34e+07    21000000       NA         NA         NA
## MSQ4_AG_152             4.49e+06     3160000       NA         NA         NA
## MSQ4_AG_153             3.20e+07     3880000       NA         NA         NA
## MSQ4_AG_154             1.88e+07     4070000       NA         NA         NA
## MSQ4_AG_155                   NA          NA       NA         NA         NA
## MSQ4_AG_156                   NA          NA       NA         NA         NA
## MSQ4_AG_157                   NA          NA       NA         NA         NA
## MSQ4_AG_159                   NA          NA       NA         NA         NA
## MSQ4_AG_160             6.56e+06      110000       NA         NA         NA
## MSQ4_AG_161             3.85e+07      437000       NA         NA         NA
## MSQ4_AG_162             5.46e+07      302000       NA         NA         NA
## MSQ4_AG_163             2.98e+07      132000       NA         NA         NA
## MSQ4_AG_164                   NA          NA       NA         NA         NA
## MSQ4_AG_165                   NA          NA       NA         NA         NA
## MSQ4_AG_167                   NA          NA       NA         NA         NA
## MSQ4_AG_168                   NA          NA       NA         NA         NA
## MSQ4_AG_169             2.29e+07     3220000       NA         NA         NA
## MSQ4_AG_170             2.07e+07     1470000       NA         NA         NA
## MSQ4_AG_171             2.00e+07     3530000       NA         NA         NA
## MSQ4_AG_172             2.55e+07     5230000       NA         NA         NA
## MSQ4_AG_173             1.24e+07     7940000       NA         NA         NA
## MSQ4_AG_175             4.82e+07    10200000       NA         NA         NA
## MSQ4_AG_176             2.43e+07     1130000       NA         NA         NA
## MSQ4_AG_177             3.24e+07     1190000       NA         NA         NA
## MSQ4_AG_178             1.44e+07      315000       NA         NA         NA
## MSQ4_AG_179             2.60e+07     3570000       NA         NA         NA
## MSQ4_AG_180             1.90e+07     1010000       NA         NA         NA
## MSQ4_AG_181             4.23e+07      142000       NA         NA         NA
## MSQ4_AG_183             7.17e+07      661000       NA         NA         NA
## MSQ4_AG_184             2.71e+07          NA       NA         NA         NA
## MSQ4_AG_185                   NA          NA       NA         NA         NA
## MSQ4_AG_186                   NA          NA       NA         NA         NA
## MSQ4_AG_187                   NA          NA       NA         NA         NA
## MSQ4_AG_188                   NA          NA       NA         NA         NA
## MSQ4_AG_189                   NA          NA       NA         NA         NA
## MSQ4_AG_191                   NA          NA       NA         NA         NA
## MSQ4_AG_192                   NA          NA       NA         NA         NA
## MSQ4_AG_193                   NA          NA       NA         NA         NA
## MSQ4_AG_194                   NA          NA       NA         NA         NA
## MSQ4_AG_195                   NA          NA       NA         NA         NA
## MSQ4_AG_196                   NA          NA       NA         NA         NA
## MSQ4_AG_197                   NA          NA       NA         NA         NA
## MSQ4_AG_198                   NA          NA       NA         NA         NA
## MSQ4_AG_201                   NA          NA       NA         NA         NA
## MSQ4_AG_205                   NA          NA       NA         NA         NA
## MSQ4_AG_207                   NA          NA       NA         NA         NA
## MSQ4_AG_208                   NA          NA       NA         NA         NA
## MSQ4_AG_209                   NA          NA       NA         NA         NA
## MSQ4_AG_210                   NA          NA       NA         NA         NA
## MSQ4_AG_211                   NA          NA       NA         NA         NA
## MSQ4_AG_212                   NA          NA       NA         NA         NA
## MSQ4_AG_213                   NA          NA       NA         NA         NA
## MSQ4_AG_215                   NA          NA       NA         NA         NA
## MSQ4_AG_216                   NA          NA       NA         NA         NA
## MSQ4_AG_217                   NA          NA       NA         NA         NA
## MSQ4_AG_218                   NA          NA       NA         NA         NA
## MSQ4_AG_219                   NA          NA       NA         NA         NA
## MSQ4_AG_220                   NA          NA       NA         NA         NA
## MSQ4_AG_221                   NA          NA       NA         NA         NA
## MSQ4_AG_223                   NA          NA       NA         NA         NA
## MSQ4_AG_224                   NA          NA       NA         NA         NA
## MSQ4_AG_225                   NA          NA       NA         NA         NA
## MSQ4_AG_226                   NA          NA       NA         NA         NA
## MSQ4_AG_227                   NA          NA       NA         NA         NA
## MSQ4_AG_228                   NA          NA       NA         NA         NA
## MSQ4_AG_229                   NA          NA       NA         NA         NA
## MSQ4_AG_230                   NA          NA       NA         NA         NA
## MSQ4_AG_231                   NA          NA       NA         NA         NA
## MSQ4_AG_232                   NA          NA       NA         NA         NA
## MSQ4_AG_233                   NA          NA       NA         NA         NA
## MSQ4_AG_234                   NA          NA       NA         NA         NA
## MSQ4_AG_235                   NA          NA       NA         NA         NA
## MSQ4_AG_236                   NA          NA       NA         NA         NA
## MSQ4_AG_237                   NA          NA       NA         NA         NA
## MSQ4_AG_238                   NA          NA       NA         NA         NA
## MSQ4_AG_239                   NA          NA       NA         NA         NA
## MSQ4_AG_240                   NA          NA       NA         NA         NA
## MSQ4_AG_241                   NA          NA       NA         NA         NA
## MSQ4_AG_242                   NA          NA       NA         NA         NA
## MSQ4_AG_243                   NA          NA       NA         NA         NA
## MSQ4_AG_244                   NA          NA       NA         NA         NA
## MSQ4_AG_245                   NA          NA       NA         NA         NA
## MSQ4_AG_247                   NA          NA       NA         NA         NA
## MSQ4_AG_248                   NA          NA       NA         NA         NA
## MSQ4_AG_249                   NA          NA       NA         NA         NA
## MSQ6_AG_250                   NA          NA       NA         NA         NA
## MSQ6_AG_251                   NA          NA       NA         NA         NA
## MSQ6_AG_252                   NA          NA       NA         NA         NA
## MSQ6_AG_253                   NA          NA       NA         NA         NA
## MSQ6_AG_254                   NA          NA       NA         NA         NA
## MSQ6_AG_255                   NA          NA       NA         NA         NA
## MSQ6_AG_256                   NA          NA       NA         NA         NA
## MSQ6_AG_257                   NA          NA       NA         NA         NA
## MSQ6_AG_258                   NA          NA       NA         NA         NA
## MSQ6_AG_259                   NA          NA       NA         NA         NA
## MSQ6_AG_260                   NA          NA       NA         NA         NA
## MSQ6_AG_261                   NA          NA       NA         NA         NA
## MSQ6_AG_262                   NA          NA       NA         NA         NA
## MSQ6_AG_263                   NA          NA       NA         NA         NA
## MSQ6_AG_264                   NA          NA       NA         NA         NA
## MSQ6_AG_265                   NA          NA       NA         NA         NA
## MSQ6_AG_266                   NA          NA       NA         NA         NA
## MSQ6_AG_267                   NA          NA       NA         NA         NA
## MSQ6_AG_268                   NA          NA       NA         NA         NA
## MSQ6_AG_269                   NA          NA       NA         NA         NA
## MSQ6_AG_270                   NA          NA       NA         NA         NA
## MSQ6_AG_271                   NA          NA       NA         NA         NA
## MSQ6_AG_272                   NA          NA       NA         NA         NA
## MSQ6_AG_273                   NA          NA       NA         NA         NA
## MSQ6_AG_274                   NA          NA       NA         NA         NA
## MSQ6_AG_275                   NA          NA       NA         NA         NA
## MSQ6_AG_276                   NA          NA       NA         NA         NA
## MSQ6_AG_277                   NA          NA       NA         NA         NA
## MSQ6_AG_278                   NA          NA       NA         NA         NA
## MSQ6_AG_279                   NA          NA       NA         NA         NA
## MSQ6_AG_281                   NA          NA       NA         NA         NA
## MSQ6_AG_282                   NA          NA       NA         NA         NA
## MSQ6_AG_283                   NA          NA       NA         NA         NA
## MSQ6_AG_284                   NA          NA       NA         NA         NA
## MSQ6_AG_285                   NA          NA       NA         NA         NA
## MSQ6_AG_286                   NA          NA       NA         NA         NA
## MSQ6_AG_287                   NA          NA       NA         NA         NA
## MSQ6_AG_288                   NA          NA       NA         NA         NA
## MSQ6_AG_290                   NA          NA       NA         NA         NA
## MSQ6_AG_291                   NA          NA       NA         NA         NA
## MSQ6_AG_292                   NA          NA       NA         NA         NA
## MSQ6_AG_293                   NA          NA       NA         NA         NA
## MSQ6_AG_294                   NA          NA       NA         NA         NA
## MSQ6_AG_295                   NA          NA       NA         NA         NA
## MSQ6_AG_296                   NA          NA       NA         NA         NA
## MSQ6_AG_298                   NA          NA       NA         NA         NA
## MSQ6_AG_299                   NA          NA       NA         NA         NA
## MSQ6_AG_300                   NA          NA       NA         NA         NA
## MSQ6_AG_301                   NA          NA       NA         NA         NA
## MSQ6_AG_302                   NA          NA       NA         NA         NA
## MSQ6_AG_303                   NA          NA       NA         NA         NA
## MSQ6_AG_304                   NA          NA       NA         NA         NA
## MSQ6_AG_305                   NA          NA       NA         NA         NA
## MSQ6_AG_306                   NA          NA       NA         NA         NA
## MSQ6_AG_307                   NA          NA       NA         NA         NA
## MSQ6_AG_308                   NA          NA       NA         NA         NA
## MSQ6_AG_309                   NA          NA       NA         NA         NA
## MSQ6_AG_310                   NA          NA       NA         NA         NA
## MSQ6_AG_311                   NA          NA       NA         NA         NA
## MSQ6_AG_312                   NA          NA       NA         NA         NA
## MSQ6_AG_313                   NA          NA       NA         NA         NA
## MSQ6_AG_314                   NA          NA       NA         NA         NA
## MSQ6_AG_315                   NA          NA       NA         NA         NA
## MSQ6_AG_316                   NA          NA       NA         NA         NA
## MSQ6_AG_317                   NA          NA       NA         NA         NA
## MSQ6_AG_318                   NA          NA       NA         NA         NA
## MSQ6_AG_319                   NA          NA       NA         NA         NA
## MSQ6_AG_320                   NA          NA       NA         NA         NA
## MSQ6_AG_321                   NA          NA       NA         NA         NA
## MSQ6_AG_322                   NA          NA       NA         NA         NA
## MSQ6_AG_323                   NA          NA       NA         NA         NA
## MSQ6_AG_324                   NA          NA       NA         NA         NA
## MSQ6_AG_325                   NA          NA       NA         NA         NA
## MSQ6_AG_326                   NA          NA       NA         NA         NA
## MSQ6_AG_327                   NA          NA       NA         NA         NA
## MSQ6_AG_328                   NA          NA       NA         NA         NA
## MSQ6_AG_329                   NA          NA       NA         NA         NA
## MSQ6_AG_330                   NA          NA       NA         NA         NA
## MSQ6_AG_331                   NA          NA       NA         NA         NA
## MSQ6_AG_332                   NA          NA       NA         NA         NA
## MSQ6_AG_333                   NA          NA       NA         NA         NA
## MSQ6_AG_334                   NA          NA       NA         NA         NA
## MSQ6_AG_335                   NA          NA       NA         NA         NA
## MSQ6_AG_336                   NA          NA       NA         NA         NA
## MSQ6_AG_337                   NA          NA       NA         NA         NA
## MSQ6_AG_339                   NA          NA       NA         NA         NA
## MSQ6_AG_340                   NA          NA       NA         NA         NA
## MSQ6_AG_341                   NA          NA       NA         NA         NA
## MSQ6_AG_342                   NA          NA       NA         NA         NA
## MSQ6_AG_343                   NA          NA       NA         NA         NA
## MSQ6_AG_344                   NA          NA       NA         NA         NA
## MSQ6_AG_346                   NA          NA       NA         NA         NA
## MSQ6_AG_347                   NA          NA       NA         NA         NA
## MSQ6_AG_348                   NA          NA       NA         NA         NA
## MSQ6_AG_349                   NA          NA       NA         NA         NA
## MSQ6_AG_350                   NA          NA       NA         NA         NA
## MSQ6_AG_351                   NA          NA       NA         NA         NA
## MSQ6_AG_352                   NA          NA       NA         NA         NA
## MSQ6_AG_353                   NA          NA       NA         NA         NA
## MSQ6_AG_354                   NA          NA       NA         NA         NA
## MSQ6_AG_355                   NA          NA       NA         NA         NA
## MSQ6_AG_356                   NA          NA       NA         NA         NA
## MSQ6_AG_357                   NA          NA       NA         NA         NA
## MSQ6_AG_358                   NA          NA       NA         NA         NA
## MSQ6_AG_359                   NA          NA       NA         NA         NA
## MSQ6_AG_360                   NA          NA       NA         NA         NA
## MSQ6_AG_361                   NA          NA       NA         NA         NA
## MSQ6_AG_362                   NA          NA       NA         NA         NA
## MSQ6_AG_363                   NA          NA       NA         NA         NA
## MSQ6_AG_364                   NA          NA       NA         NA         NA
## MSQ6_AG_365                   NA          NA       NA         NA         NA
## MSQ6_AG_366                   NA          NA       NA         NA         NA
## MSQ6_AG_367                   NA          NA       NA         NA         NA
## MSQ6_AG_368                   NA          NA       NA         NA         NA
## MSQ6_AG_369                   NA          NA       NA         NA         NA
## MSQ6_AG_370                   NA          NA       NA         NA         NA
## MSQ6_AG_371                   NA          NA       NA         NA         NA
## MSQ6_AG_372                   NA          NA       NA         NA         NA
## MSQ6_AG_373                   NA          NA       NA         NA         NA
## MSQ6_AG_374                   NA          NA       NA         NA         NA
## MSQ6_AG_375                   NA          NA       NA         NA         NA
## MSQ6_AG_376                   NA          NA       NA         NA         NA
## MSQ6_AG_377                   NA          NA       NA         NA         NA
## MSQ6_AG_378                   NA          NA       NA         NA         NA
## MSQ6_AG_379                   NA          NA       NA         NA         NA
## MSQ6_AG_380                   NA          NA       NA         NA         NA
## MSQ6_AG_381                   NA          NA       NA         NA         NA
## MSQ6_AG_382                   NA          NA       NA         NA         NA
## MSQ6_AG_383                   NA          NA       NA         NA         NA
## MSQ6_AG_384                   NA          NA       NA         NA         NA
## MSQ6_AG_385                   NA          NA       NA         NA         NA
## MSQ6_AG_386                   NA          NA       NA         NA         NA
## MSQ6_AG_387                   NA          NA       NA         NA         NA
## MSQ6_AG_388             7.36e+06       38400       NA         NA         NA
## MSQ6_AG_389             1.84e+07      127000       NA         NA         NA
## MSQ6_AG_390             1.69e+07      674000       NA         NA         NA
## MSQ6_AG_391             1.53e+06       21500       NA         NA         NA
## MSQ6_AG_392             2.47e+07      483000       NA         NA         NA
## MSQ6_AG_393             3.92e+07      249000       NA         NA         NA
## MSQ6_AG_394             3.91e+07      252000       NA         NA         NA
## MSQ6_AG_395             1.84e+07      208000       NA         NA         NA
## MSQ6_AG_396             1.43e+07      201000       NA         NA         NA
## MSQ6_AG_397             1.26e+07       90100       NA         NA         NA
## MSQ6_AG_398             4.46e+07    18100000       NA         NA         NA
## MSQ6_AG_399             3.79e+07      483000       NA         NA         NA
## MSQ6_AG_400             2.40e+07      147000       NA         NA         NA
## MSQ6_AG_401             3.56e+07      968000       NA         NA         NA
## MSQ6_AG_402             1.43e+07      725000       NA         NA         NA
## MSQ6_AG_403             5.93e+07     9750000       NA         NA         NA
## MSQ6_AG_404             1.00e+08     1070000       NA         NA         NA
## MSQ6_AG_405             3.52e+07      221000       NA         NA         NA
## MSQ6_AG_406             2.61e+07      255000       NA         NA         NA
## MSQ6_AG_407             2.44e+07      268000       NA         NA         NA
## MSQ6_AG_408             3.81e+07      924000       NA         NA         NA
## MSQ6_AG_409             4.68e+07      365000       NA         NA         NA
## MSQ6_AG_410             1.06e+08      617000       NA         NA         NA
## MSQ6_AG_411             9.34e+07      602000       NA         NA         NA
## MSQ6_AG_412             2.85e+07     6850000       NA         NA         NA
## MSQ6_AG_413             4.41e+07      333000       NA         NA         NA
## MSQ6_AG_414             5.80e+07      556000       NA         NA         NA
## MSQ6_AG_415             5.76e+07      231000       NA         NA         NA
## MSQ6_AG_416             3.03e+06      192000       NA         NA         NA
## MSQ6_AG_417             6.15e+07      497000       NA         NA         NA
## MSQ6_AG_418             5.52e+07     4270000       NA         NA         NA
## MSQ6_AG_419             2.93e+07      402000       NA         NA         NA
## MSQ6_AG_420             4.44e+07      592000       NA         NA         NA
## MSQ6_AG_421             9.94e+07      535000       NA         NA         NA
## MSQ6_AG_422             2.62e+07     1810000       NA         NA         NA
## MSQ6_AG_423             8.83e+07      852000       NA         NA         NA
## MSQ6_AG_424             9.22e+07      240000       NA         NA         NA
## MSQ6_AG_425             6.42e+07      557000       NA         NA         NA
## MSQ6_AG_426             9.21e+07     2430000       NA         NA         NA
## MSQ6_AG_427             1.11e+08    41900000       NA         NA         NA
## MSQ6_AG_428             9.63e+07    15000000       NA         NA         NA
## MSQ6_AG_429             2.97e+08    18000000       NA         NA         NA
## MSQ6_AG_430             1.20e+08    17400000       NA         NA         NA
## MSQ6_AG_431             5.07e+08    83300000       NA         NA         NA
## MSQ6_AG_432             1.04e+08    12500000       NA         NA         NA
## MSQ6_AG_433             3.08e+07     5340000       NA         NA         NA
## MSQ6_AG_434             1.86e+08    31400000       NA         NA         NA
## MSQ6_AG_490                   NA          NA       NA         NA         NA
## MSQ6_AG_491                   NA          NA       NA         NA         NA
## MSQ6_AG_492                   NA          NA       NA         NA         NA
## MSQ6_AG_493                   NA          NA       NA         NA         NA
## MSQ6_AG_494                   NA          NA       NA         NA         NA
## MSQ6_AG_495                   NA          NA       NA         NA         NA
## MSQ6_AG_496                   NA          NA       NA         NA         NA
## MSQ6_AG_497                   NA          NA       NA         NA         NA
## MSQ6_AG_498                   NA          NA       NA         NA         NA
## MSQ6_AG_499                   NA          NA       NA         NA         NA
## MSQ6_AG_500                   NA          NA       NA         NA         NA
## MSQ6_AG_501                   NA          NA       NA         NA         NA
## MSQ6_AG_502                   NA          NA       NA         NA         NA
## MSQ6_AG_503                   NA          NA       NA         NA         NA
## MSQ6_AG_504                   NA          NA       NA         NA         NA
## MSQ6_AG_505                   NA          NA       NA         NA         NA
## MSQ6_AG_506                   NA          NA       NA         NA         NA
## MSQ6_AG_509                   NA          NA       NA         NA         NA
## MSQ6_AG_511                   NA          NA       NA         NA         NA
## MSQ6_AG_512                   NA          NA       NA         NA         NA
## MSQ6_AG_514                   NA          NA       NA         NA         NA
## MSQ6_AG_515                   NA          NA       NA         NA         NA
## MSQ6_AG_516                   NA          NA       NA         NA         NA
## MSQ6_AG_517                   NA          NA       NA         NA         NA
## MSQ6_AG_518                   NA          NA       NA         NA         NA
## MSQ6_AG_519                   NA          NA       NA         NA         NA
## MSQ6_AG_520                   NA          NA       NA         NA         NA
## MSQ6_AG_521                   NA          NA       NA         NA         NA
## MSQ6_AG_522                   NA          NA       NA         NA         NA
## MSQ6_AG_523                   NA          NA       NA         NA         NA
## MSQ6_AG_524                   NA          NA       NA         NA         NA
## MSQ6_AG_525                   NA          NA       NA         NA         NA
## MSQ6_AG_526                   NA          NA       NA         NA         NA
## MSQ6_AG_527                   NA          NA       NA         NA         NA
## MSQ6_AG_528                   NA          NA       NA         NA         NA
## MSQ6_AG_529                   NA          NA       NA         NA         NA
## MSQ6_AG_530                   NA          NA       NA         NA         NA
##             Spermidine Spermine B24_Acetate B24_Propionate B24_Butyrate
## MSQ4_AG_012         NA       NA          NA             NA           NA
## MSQ4_AG_013         NA       NA          NA             NA           NA
## MSQ4_AG_014         NA       NA          NA             NA           NA
## MSQ4_AG_015         NA       NA          NA             NA           NA
## MSQ4_AG_016         NA       NA          NA             NA           NA
## MSQ4_AG_017         NA       NA          NA             NA           NA
## MSQ4_AG_018         NA       NA          NA             NA           NA
## MSQ4_AG_019         NA       NA          NA             NA           NA
## MSQ4_AG_020         NA       NA          NA             NA           NA
## MSQ4_AG_021         NA       NA          NA             NA           NA
## MSQ4_AG_022         NA       NA          NA             NA           NA
## MSQ4_AG_023         NA       NA          NA             NA           NA
## MSQ4_AG_024         NA       NA          NA             NA           NA
## MSQ4_AG_025         NA       NA          NA             NA           NA
## MSQ4_AG_027         NA       NA          NA             NA           NA
## MSQ4_AG_028         NA       NA          NA             NA           NA
## MSQ4_AG_029         NA       NA          NA             NA           NA
## MSQ4_AG_030         NA       NA          NA             NA           NA
## MSQ4_AG_031         NA       NA          NA             NA           NA
## MSQ4_AG_033         NA       NA          NA             NA           NA
## MSQ4_AG_034         NA       NA          NA             NA           NA
## MSQ4_AG_035         NA       NA          NA             NA           NA
## MSQ4_AG_036         NA       NA          NA             NA           NA
## MSQ4_AG_037         NA       NA          NA             NA           NA
## MSQ4_AG_038         NA       NA          NA             NA           NA
## MSQ4_AG_039         NA       NA          NA             NA           NA
## MSQ4_AG_040         NA       NA          NA             NA           NA
## MSQ4_AG_041         NA       NA          NA             NA           NA
## MSQ4_AG_042         NA       NA          NA             NA           NA
## MSQ4_AG_043         NA       NA          NA             NA           NA
## MSQ4_AG_044         NA       NA          NA             NA           NA
## MSQ4_AG_045         NA       NA          NA             NA           NA
## MSQ4_AG_046         NA       NA          NA             NA           NA
## MSQ4_AG_047         NA       NA          NA             NA           NA
## MSQ4_AG_049         NA       NA          NA             NA           NA
## MSQ4_AG_050         NA       NA          NA             NA           NA
## MSQ4_AG_051         NA       NA          NA             NA           NA
## MSQ4_AG_052         NA       NA          NA             NA           NA
## MSQ4_AG_053         NA       NA          NA             NA           NA
## MSQ4_AG_054         NA       NA          NA             NA           NA
## MSQ4_AG_055         NA       NA          NA             NA           NA
## MSQ4_AG_056         NA       NA          NA             NA           NA
## MSQ4_AG_057         NA       NA          NA             NA           NA
## MSQ4_AG_058         NA       NA          NA             NA           NA
## MSQ4_AG_059         NA       NA          NA             NA           NA
## MSQ4_AG_062         NA       NA          NA             NA           NA
## MSQ4_AG_064         NA       NA          NA             NA           NA
## MSQ4_AG_065         NA       NA          NA             NA           NA
## MSQ4_AG_072         NA       NA          NA             NA           NA
## MSQ4_AG_076         NA       NA          NA             NA           NA
## MSQ4_AG_077         NA       NA          NA             NA           NA
## MSQ4_AG_080         NA       NA          NA             NA           NA
## MSQ4_AG_083         NA       NA          NA             NA           NA
## MSQ4_AG_084         NA       NA          NA             NA           NA
## MSQ4_AG_085         NA       NA          NA             NA           NA
## MSQ4_AG_086         NA       NA          NA             NA           NA
## MSQ4_AG_087         NA       NA          NA             NA           NA
## MSQ4_AG_088         NA       NA          NA             NA           NA
## MSQ4_AG_089         NA       NA          NA             NA           NA
## MSQ4_AG_091         NA       NA          NA             NA           NA
## MSQ4_AG_092         NA       NA          NA             NA           NA
## MSQ4_AG_093         NA       NA          NA             NA           NA
## MSQ4_AG_095         NA       NA          NA             NA           NA
## MSQ4_AG_096         NA       NA          NA             NA           NA
## MSQ4_AG_101         NA       NA          NA             NA           NA
## MSQ4_AG_103         NA       NA          NA             NA           NA
## MSQ4_AG_104         NA       NA          NA             NA           NA
## MSQ4_AG_105         NA       NA          NA             NA           NA
## MSQ4_AG_106         NA       NA          NA             NA           NA
## MSQ4_AG_107         NA       NA          NA             NA           NA
## MSQ4_AG_108         NA       NA          NA             NA           NA
## MSQ4_AG_109         NA       NA          NA             NA           NA
## MSQ4_AG_111         NA       NA          NA             NA           NA
## MSQ4_AG_112         NA       NA          NA             NA           NA
## MSQ4_AG_113         NA       NA          NA             NA           NA
## MSQ4_AG_117         NA       NA          NA             NA           NA
## MSQ4_AG_119         NA       NA          NA             NA           NA
## MSQ4_AG_121         NA       NA          NA             NA           NA
## MSQ4_AG_125         NA       NA          NA             NA           NA
## MSQ4_AG_127         NA       NA          NA             NA           NA
## MSQ4_AG_128         NA       NA          NA             NA           NA
## MSQ4_AG_129         NA       NA          NA             NA           NA
## MSQ4_AG_130         NA       NA          NA             NA           NA
## MSQ4_AG_131         NA       NA          NA             NA           NA
## MSQ4_AG_132         NA       NA          NA             NA           NA
## MSQ4_AG_133         NA       NA          NA             NA           NA
## MSQ4_AG_135         NA       NA          NA             NA           NA
## MSQ4_AG_136         NA       NA          NA             NA           NA
## MSQ4_AG_141         NA       NA          NA             NA           NA
## MSQ4_AG_143         NA       NA          NA             NA           NA
## MSQ4_AG_144         NA       NA          NA             NA           NA
## MSQ4_AG_145         NA       NA          NA             NA           NA
## MSQ4_AG_148         NA       NA          NA             NA           NA
## MSQ4_AG_149         NA       NA          NA             NA           NA
## MSQ4_AG_151         NA       NA          NA             NA           NA
## MSQ4_AG_152         NA       NA          NA             NA           NA
## MSQ4_AG_153         NA       NA          NA             NA           NA
## MSQ4_AG_154         NA       NA          NA             NA           NA
## MSQ4_AG_155         NA       NA          NA             NA           NA
## MSQ4_AG_156         NA       NA          NA             NA           NA
## MSQ4_AG_157         NA       NA          NA             NA           NA
## MSQ4_AG_159         NA       NA          NA             NA           NA
## MSQ4_AG_160         NA       NA          NA             NA           NA
## MSQ4_AG_161         NA       NA          NA             NA           NA
## MSQ4_AG_162         NA       NA          NA             NA           NA
## MSQ4_AG_163         NA       NA          NA             NA           NA
## MSQ4_AG_164         NA       NA          NA             NA           NA
## MSQ4_AG_165         NA       NA          NA             NA           NA
## MSQ4_AG_167         NA       NA          NA             NA           NA
## MSQ4_AG_168         NA       NA          NA             NA           NA
## MSQ4_AG_169         NA       NA          NA             NA           NA
## MSQ4_AG_170         NA       NA          NA             NA           NA
## MSQ4_AG_171         NA       NA          NA             NA           NA
## MSQ4_AG_172         NA       NA          NA             NA           NA
## MSQ4_AG_173         NA       NA          NA             NA           NA
## MSQ4_AG_175         NA       NA          NA             NA           NA
## MSQ4_AG_176         NA       NA          NA             NA           NA
## MSQ4_AG_177         NA       NA          NA             NA           NA
## MSQ4_AG_178         NA       NA          NA             NA           NA
## MSQ4_AG_179         NA       NA          NA             NA           NA
## MSQ4_AG_180         NA       NA          NA             NA           NA
## MSQ4_AG_181         NA       NA          NA             NA           NA
## MSQ4_AG_183         NA       NA          NA             NA           NA
## MSQ4_AG_184         NA       NA          NA             NA           NA
## MSQ4_AG_185         NA       NA          NA             NA           NA
## MSQ4_AG_186         NA       NA          NA             NA           NA
## MSQ4_AG_187         NA       NA          NA             NA           NA
## MSQ4_AG_188         NA       NA          NA             NA           NA
## MSQ4_AG_189         NA       NA          NA             NA           NA
## MSQ4_AG_191         NA       NA          NA             NA           NA
## MSQ4_AG_192         NA       NA          NA             NA           NA
## MSQ4_AG_193         NA       NA          NA             NA           NA
## MSQ4_AG_194         NA       NA          NA             NA           NA
## MSQ4_AG_195         NA       NA          NA             NA           NA
## MSQ4_AG_196         NA       NA          NA             NA           NA
## MSQ4_AG_197         NA       NA          NA             NA           NA
## MSQ4_AG_198         NA       NA          NA             NA           NA
## MSQ4_AG_201         NA       NA          NA             NA           NA
## MSQ4_AG_205         NA       NA          NA             NA           NA
## MSQ4_AG_207         NA       NA          NA             NA           NA
## MSQ4_AG_208         NA       NA          NA             NA           NA
## MSQ4_AG_209         NA       NA          NA             NA           NA
## MSQ4_AG_210         NA       NA          NA             NA           NA
## MSQ4_AG_211         NA       NA          NA             NA           NA
## MSQ4_AG_212         NA       NA          NA             NA           NA
## MSQ4_AG_213         NA       NA          NA             NA           NA
## MSQ4_AG_215         NA       NA          NA             NA           NA
## MSQ4_AG_216         NA       NA          NA             NA           NA
## MSQ4_AG_217         NA       NA          NA             NA           NA
## MSQ4_AG_218         NA       NA          NA             NA           NA
## MSQ4_AG_219         NA       NA          NA             NA           NA
## MSQ4_AG_220         NA       NA          NA             NA           NA
## MSQ4_AG_221         NA       NA          NA             NA           NA
## MSQ4_AG_223         NA       NA          NA             NA           NA
## MSQ4_AG_224         NA       NA          NA             NA           NA
## MSQ4_AG_225         NA       NA          NA             NA           NA
## MSQ4_AG_226         NA       NA          NA             NA           NA
## MSQ4_AG_227         NA       NA          NA             NA           NA
## MSQ4_AG_228         NA       NA          NA             NA           NA
## MSQ4_AG_229         NA       NA          NA             NA           NA
## MSQ4_AG_230         NA       NA          NA             NA           NA
## MSQ4_AG_231         NA       NA          NA             NA           NA
## MSQ4_AG_232         NA       NA          NA             NA           NA
## MSQ4_AG_233         NA       NA          NA             NA           NA
## MSQ4_AG_234         NA       NA          NA             NA           NA
## MSQ4_AG_235         NA       NA          NA             NA           NA
## MSQ4_AG_236         NA       NA          NA             NA           NA
## MSQ4_AG_237         NA       NA          NA             NA           NA
## MSQ4_AG_238         NA       NA          NA             NA           NA
## MSQ4_AG_239         NA       NA          NA             NA           NA
## MSQ4_AG_240         NA       NA          NA             NA           NA
## MSQ4_AG_241         NA       NA          NA             NA           NA
## MSQ4_AG_242         NA       NA          NA             NA           NA
## MSQ4_AG_243         NA       NA          NA             NA           NA
## MSQ4_AG_244         NA       NA          NA             NA           NA
## MSQ4_AG_245         NA       NA          NA             NA           NA
## MSQ4_AG_247         NA       NA          NA             NA           NA
## MSQ4_AG_248         NA       NA          NA             NA           NA
## MSQ4_AG_249         NA       NA          NA             NA           NA
## MSQ6_AG_250         NA       NA          NA             NA           NA
## MSQ6_AG_251         NA       NA          NA             NA           NA
## MSQ6_AG_252         NA       NA          NA             NA           NA
## MSQ6_AG_253         NA       NA          NA             NA           NA
## MSQ6_AG_254         NA       NA          NA             NA           NA
## MSQ6_AG_255         NA       NA          NA             NA           NA
## MSQ6_AG_256         NA       NA          NA             NA           NA
## MSQ6_AG_257         NA       NA          NA             NA           NA
## MSQ6_AG_258         NA       NA          NA             NA           NA
## MSQ6_AG_259         NA       NA          NA             NA           NA
## MSQ6_AG_260         NA       NA          NA             NA           NA
## MSQ6_AG_261         NA       NA          NA             NA           NA
## MSQ6_AG_262         NA       NA          NA             NA           NA
## MSQ6_AG_263         NA       NA          NA             NA           NA
## MSQ6_AG_264         NA       NA          NA             NA           NA
## MSQ6_AG_265         NA       NA          NA             NA           NA
## MSQ6_AG_266         NA       NA          NA             NA           NA
## MSQ6_AG_267         NA       NA          NA             NA           NA
## MSQ6_AG_268         NA       NA          NA             NA           NA
## MSQ6_AG_269         NA       NA          NA             NA           NA
## MSQ6_AG_270         NA       NA          NA             NA           NA
## MSQ6_AG_271         NA       NA          NA             NA           NA
## MSQ6_AG_272         NA       NA          NA             NA           NA
## MSQ6_AG_273         NA       NA          NA             NA           NA
## MSQ6_AG_274         NA       NA          NA             NA           NA
## MSQ6_AG_275         NA       NA          NA             NA           NA
## MSQ6_AG_276         NA       NA          NA             NA           NA
## MSQ6_AG_277         NA       NA          NA             NA           NA
## MSQ6_AG_278         NA       NA          NA             NA           NA
## MSQ6_AG_279         NA       NA          NA             NA           NA
## MSQ6_AG_281         NA       NA          NA             NA           NA
## MSQ6_AG_282         NA       NA          NA             NA           NA
## MSQ6_AG_283         NA       NA          NA             NA           NA
## MSQ6_AG_284         NA       NA          NA             NA           NA
## MSQ6_AG_285         NA       NA          NA             NA           NA
## MSQ6_AG_286         NA       NA          NA             NA           NA
## MSQ6_AG_287         NA       NA          NA             NA           NA
## MSQ6_AG_288         NA       NA          NA             NA           NA
## MSQ6_AG_290         NA       NA          NA             NA           NA
## MSQ6_AG_291         NA       NA          NA             NA           NA
## MSQ6_AG_292         NA       NA          NA             NA           NA
## MSQ6_AG_293         NA       NA          NA             NA           NA
## MSQ6_AG_294         NA       NA          NA             NA           NA
## MSQ6_AG_295         NA       NA          NA             NA           NA
## MSQ6_AG_296         NA       NA          NA             NA           NA
## MSQ6_AG_298         NA       NA          NA             NA           NA
## MSQ6_AG_299         NA       NA          NA             NA           NA
## MSQ6_AG_300         NA       NA          NA             NA           NA
## MSQ6_AG_301         NA       NA          NA             NA           NA
## MSQ6_AG_302         NA       NA          NA             NA           NA
## MSQ6_AG_303         NA       NA          NA             NA           NA
## MSQ6_AG_304         NA       NA          NA             NA           NA
## MSQ6_AG_305         NA       NA          NA             NA           NA
## MSQ6_AG_306         NA       NA          NA             NA           NA
## MSQ6_AG_307         NA       NA          NA             NA           NA
## MSQ6_AG_308         NA       NA          NA             NA           NA
## MSQ6_AG_309         NA       NA          NA             NA           NA
## MSQ6_AG_310         NA       NA          NA             NA           NA
## MSQ6_AG_311         NA       NA          NA             NA           NA
## MSQ6_AG_312         NA       NA          NA             NA           NA
## MSQ6_AG_313         NA       NA          NA             NA           NA
## MSQ6_AG_314         NA       NA          NA             NA           NA
## MSQ6_AG_315         NA       NA          NA             NA           NA
## MSQ6_AG_316         NA       NA          NA             NA           NA
## MSQ6_AG_317         NA       NA          NA             NA           NA
## MSQ6_AG_318         NA       NA          NA             NA           NA
## MSQ6_AG_319         NA       NA          NA             NA           NA
## MSQ6_AG_320         NA       NA          NA             NA           NA
## MSQ6_AG_321         NA       NA          NA             NA           NA
## MSQ6_AG_322         NA       NA          NA             NA           NA
## MSQ6_AG_323         NA       NA          NA             NA           NA
## MSQ6_AG_324         NA       NA          NA             NA           NA
## MSQ6_AG_325         NA       NA          NA             NA           NA
## MSQ6_AG_326         NA       NA          NA             NA           NA
## MSQ6_AG_327         NA       NA          NA             NA           NA
## MSQ6_AG_328         NA       NA          NA             NA           NA
## MSQ6_AG_329         NA       NA          NA             NA           NA
## MSQ6_AG_330         NA       NA          NA             NA           NA
## MSQ6_AG_331         NA       NA          NA             NA           NA
## MSQ6_AG_332         NA       NA          NA             NA           NA
## MSQ6_AG_333         NA       NA          NA             NA           NA
## MSQ6_AG_334         NA       NA          NA             NA           NA
## MSQ6_AG_335         NA       NA          NA             NA           NA
## MSQ6_AG_336         NA       NA          NA             NA           NA
## MSQ6_AG_337         NA       NA          NA             NA           NA
## MSQ6_AG_339         NA       NA          NA             NA           NA
## MSQ6_AG_340         NA       NA          NA             NA           NA
## MSQ6_AG_341         NA       NA          NA             NA           NA
## MSQ6_AG_342         NA       NA          NA             NA           NA
## MSQ6_AG_343         NA       NA          NA             NA           NA
## MSQ6_AG_344         NA       NA          NA             NA           NA
## MSQ6_AG_346         NA       NA          NA             NA           NA
## MSQ6_AG_347         NA       NA          NA             NA           NA
## MSQ6_AG_348         NA       NA          NA             NA           NA
## MSQ6_AG_349         NA       NA          NA             NA           NA
## MSQ6_AG_350         NA       NA          NA             NA           NA
## MSQ6_AG_351         NA       NA          NA             NA           NA
## MSQ6_AG_352         NA       NA          NA             NA           NA
## MSQ6_AG_353         NA       NA          NA             NA           NA
## MSQ6_AG_354         NA       NA          NA             NA           NA
## MSQ6_AG_355         NA       NA          NA             NA           NA
## MSQ6_AG_356         NA       NA          NA             NA           NA
## MSQ6_AG_357         NA       NA          NA             NA           NA
## MSQ6_AG_358         NA       NA          NA             NA           NA
## MSQ6_AG_359         NA       NA          NA             NA           NA
## MSQ6_AG_360         NA       NA          NA             NA           NA
## MSQ6_AG_361         NA       NA          NA             NA           NA
## MSQ6_AG_362         NA       NA          NA             NA           NA
## MSQ6_AG_363         NA       NA          NA             NA           NA
## MSQ6_AG_364         NA       NA          NA             NA           NA
## MSQ6_AG_365         NA       NA          NA             NA           NA
## MSQ6_AG_366         NA       NA          NA             NA           NA
## MSQ6_AG_367         NA       NA          NA             NA           NA
## MSQ6_AG_368         NA       NA          NA             NA           NA
## MSQ6_AG_369         NA       NA          NA             NA           NA
## MSQ6_AG_370         NA       NA          NA             NA           NA
## MSQ6_AG_371         NA       NA          NA             NA           NA
## MSQ6_AG_372         NA       NA          NA             NA           NA
## MSQ6_AG_373         NA       NA          NA             NA           NA
## MSQ6_AG_374         NA       NA          NA             NA           NA
## MSQ6_AG_375         NA       NA          NA             NA           NA
## MSQ6_AG_376         NA       NA          NA             NA           NA
## MSQ6_AG_377         NA       NA          NA             NA           NA
## MSQ6_AG_378         NA       NA          NA             NA           NA
## MSQ6_AG_379         NA       NA          NA             NA           NA
## MSQ6_AG_380         NA       NA          NA             NA           NA
## MSQ6_AG_381         NA       NA          NA             NA           NA
## MSQ6_AG_382         NA       NA          NA             NA           NA
## MSQ6_AG_383         NA       NA          NA             NA           NA
## MSQ6_AG_384         NA       NA          NA             NA           NA
## MSQ6_AG_385         NA       NA          NA             NA           NA
## MSQ6_AG_386         NA       NA          NA             NA           NA
## MSQ6_AG_387         NA       NA          NA             NA           NA
## MSQ6_AG_388         NA       NA      35.042          6.157        4.854
## MSQ6_AG_389         NA       NA      42.309          7.751        5.909
## MSQ6_AG_390         NA       NA      42.234          6.889        5.705
## MSQ6_AG_391         NA       NA      33.957          2.847        7.643
## MSQ6_AG_392         NA       NA      52.591         14.767        7.361
## MSQ6_AG_393         NA       NA      47.967         10.082        7.990
## MSQ6_AG_394         NA       NA      44.142         10.472        8.230
## MSQ6_AG_395         NA       NA      61.325         10.666        9.577
## MSQ6_AG_396         NA       NA      49.212         10.615        8.490
## MSQ6_AG_397         NA       NA      43.532          5.885       10.808
## MSQ6_AG_398         NA       NA      43.557          6.323        7.695
## MSQ6_AG_399         NA       NA      43.178          9.706       11.763
## MSQ6_AG_400         NA       NA      43.129          7.570       11.634
## MSQ6_AG_401         NA       NA      52.961         14.023        9.507
## MSQ6_AG_402         NA       NA      45.456          8.727        7.941
## MSQ6_AG_403         NA       NA      42.449          3.745        3.872
## MSQ6_AG_404         NA       NA      44.685          8.222        8.278
## MSQ6_AG_405         NA       NA      42.383          8.581        7.299
## MSQ6_AG_406         NA       NA      42.374          6.263        7.124
## MSQ6_AG_407         NA       NA      42.785          6.478        6.081
## MSQ6_AG_408         NA       NA      38.994          6.891        5.238
## MSQ6_AG_409         NA       NA      38.766          6.084        5.346
## MSQ6_AG_410         NA       NA      41.839          4.994        5.352
## MSQ6_AG_411         NA       NA      45.235          6.486        6.220
## MSQ6_AG_412         NA       NA      42.752          5.103        5.836
## MSQ6_AG_413         NA       NA      42.962          4.288        4.468
## MSQ6_AG_414         NA       NA      42.644          5.281        3.711
## MSQ6_AG_415         NA       NA      39.557          3.525        7.333
## MSQ6_AG_416         NA       NA      32.876          1.922        1.108
## MSQ6_AG_417         NA       NA      39.641          6.139        4.911
## MSQ6_AG_418         NA       NA      45.913          8.426        7.029
## MSQ6_AG_419         NA       NA      40.378          5.077        5.000
## MSQ6_AG_420         NA       NA      44.921          8.815        9.052
## MSQ6_AG_421         NA       NA      40.332          2.766        3.897
## MSQ6_AG_422         NA       NA      39.582          7.554        8.683
## MSQ6_AG_423         NA       NA      35.387          8.234       10.618
## MSQ6_AG_424         NA       NA      41.176          7.429        7.477
## MSQ6_AG_425         NA       NA      37.701          5.163        4.153
## MSQ6_AG_426         NA       NA      29.728          0.289        2.262
## MSQ6_AG_427         NA       NA      43.787          6.130       13.191
## MSQ6_AG_428         NA       NA      41.653          4.921       13.734
## MSQ6_AG_429         NA       NA      48.924          7.750       18.915
## MSQ6_AG_430         NA       NA      43.760          6.934       18.147
## MSQ6_AG_431         NA       NA      49.028          6.462       16.515
## MSQ6_AG_432         NA       NA      46.173          5.278       14.914
## MSQ6_AG_433         NA       NA      45.913          6.601       16.091
## MSQ6_AG_434         NA       NA      46.394          6.434       16.829
## MSQ6_AG_490         NA       NA          NA             NA           NA
## MSQ6_AG_491         NA       NA          NA             NA           NA
## MSQ6_AG_492         NA       NA          NA             NA           NA
## MSQ6_AG_493         NA       NA          NA             NA           NA
## MSQ6_AG_494         NA       NA          NA             NA           NA
## MSQ6_AG_495         NA       NA          NA             NA           NA
## MSQ6_AG_496         NA       NA          NA             NA           NA
## MSQ6_AG_497         NA       NA          NA             NA           NA
## MSQ6_AG_498         NA       NA          NA             NA           NA
## MSQ6_AG_499         NA       NA          NA             NA           NA
## MSQ6_AG_500         NA       NA          NA             NA           NA
## MSQ6_AG_501         NA       NA          NA             NA           NA
## MSQ6_AG_502         NA       NA          NA             NA           NA
## MSQ6_AG_503         NA       NA          NA             NA           NA
## MSQ6_AG_504         NA       NA          NA             NA           NA
## MSQ6_AG_505         NA       NA          NA             NA           NA
## MSQ6_AG_506         NA       NA          NA             NA           NA
## MSQ6_AG_509         NA       NA          NA             NA           NA
## MSQ6_AG_511         NA       NA          NA             NA           NA
## MSQ6_AG_512         NA       NA          NA             NA           NA
## MSQ6_AG_514         NA       NA          NA             NA           NA
## MSQ6_AG_515         NA       NA          NA             NA           NA
## MSQ6_AG_516         NA       NA          NA             NA           NA
## MSQ6_AG_517         NA       NA          NA             NA           NA
## MSQ6_AG_518         NA       NA          NA             NA           NA
## MSQ6_AG_519         NA       NA          NA             NA           NA
## MSQ6_AG_520         NA       NA          NA             NA           NA
## MSQ6_AG_521         NA       NA          NA             NA           NA
## MSQ6_AG_522         NA       NA          NA             NA           NA
## MSQ6_AG_523         NA       NA          NA             NA           NA
## MSQ6_AG_524         NA       NA          NA             NA           NA
## MSQ6_AG_525         NA       NA          NA             NA           NA
## MSQ6_AG_526         NA       NA          NA             NA           NA
## MSQ6_AG_527         NA       NA          NA             NA           NA
## MSQ6_AG_528         NA       NA          NA             NA           NA
## MSQ6_AG_529         NA       NA          NA             NA           NA
## MSQ6_AG_530         NA       NA          NA             NA           NA
##             B24_Formate B24_Succinate B24_Lactate B24_BCFA succinate lactate
## MSQ4_AG_012          NA            NA          NA       NA        NA      NA
## MSQ4_AG_013          NA            NA          NA       NA        NA      NA
## MSQ4_AG_014          NA            NA          NA       NA        NA      NA
## MSQ4_AG_015          NA            NA          NA       NA        NA      NA
## MSQ4_AG_016          NA            NA          NA       NA        NA      NA
## MSQ4_AG_017          NA            NA          NA       NA        NA      NA
## MSQ4_AG_018          NA            NA          NA       NA        NA      NA
## MSQ4_AG_019          NA            NA          NA       NA        NA      NA
## MSQ4_AG_020          NA            NA          NA       NA        NA      NA
## MSQ4_AG_021          NA            NA          NA       NA        NA      NA
## MSQ4_AG_022          NA            NA          NA       NA        NA      NA
## MSQ4_AG_023          NA            NA          NA       NA        NA      NA
## MSQ4_AG_024          NA            NA          NA       NA        NA      NA
## MSQ4_AG_025          NA            NA          NA       NA        NA      NA
## MSQ4_AG_027          NA            NA          NA       NA        NA      NA
## MSQ4_AG_028          NA            NA          NA       NA        NA      NA
## MSQ4_AG_029          NA            NA          NA       NA        NA      NA
## MSQ4_AG_030          NA            NA          NA       NA        NA      NA
## MSQ4_AG_031          NA            NA          NA       NA        NA      NA
## MSQ4_AG_033          NA            NA          NA       NA        NA      NA
## MSQ4_AG_034          NA            NA          NA       NA        NA      NA
## MSQ4_AG_035          NA            NA          NA       NA        NA      NA
## MSQ4_AG_036          NA            NA          NA       NA        NA      NA
## MSQ4_AG_037          NA            NA          NA       NA        NA      NA
## MSQ4_AG_038          NA            NA          NA       NA        NA      NA
## MSQ4_AG_039          NA            NA          NA       NA        NA      NA
## MSQ4_AG_040          NA            NA          NA       NA        NA      NA
## MSQ4_AG_041          NA            NA          NA       NA        NA      NA
## MSQ4_AG_042          NA            NA          NA       NA        NA      NA
## MSQ4_AG_043          NA            NA          NA       NA        NA      NA
## MSQ4_AG_044          NA            NA          NA       NA        NA      NA
## MSQ4_AG_045          NA            NA          NA       NA        NA      NA
## MSQ4_AG_046          NA            NA          NA       NA        NA      NA
## MSQ4_AG_047          NA            NA          NA       NA        NA      NA
## MSQ4_AG_049          NA            NA          NA       NA        NA      NA
## MSQ4_AG_050          NA            NA          NA       NA        NA      NA
## MSQ4_AG_051          NA            NA          NA       NA        NA      NA
## MSQ4_AG_052          NA            NA          NA       NA        NA      NA
## MSQ4_AG_053          NA            NA          NA       NA        NA      NA
## MSQ4_AG_054          NA            NA          NA       NA        NA      NA
## MSQ4_AG_055          NA            NA          NA       NA        NA      NA
## MSQ4_AG_056          NA            NA          NA       NA        NA      NA
## MSQ4_AG_057          NA            NA          NA       NA        NA      NA
## MSQ4_AG_058          NA            NA          NA       NA        NA      NA
## MSQ4_AG_059          NA            NA          NA       NA        NA      NA
## MSQ4_AG_062          NA            NA          NA       NA        NA      NA
## MSQ4_AG_064          NA            NA          NA       NA        NA      NA
## MSQ4_AG_065          NA            NA          NA       NA        NA      NA
## MSQ4_AG_072          NA            NA          NA       NA        NA      NA
## MSQ4_AG_076          NA            NA          NA       NA        NA      NA
## MSQ4_AG_077          NA            NA          NA       NA        NA      NA
## MSQ4_AG_080          NA            NA          NA       NA        NA      NA
## MSQ4_AG_083          NA            NA          NA       NA        NA      NA
## MSQ4_AG_084          NA            NA          NA       NA        NA      NA
## MSQ4_AG_085          NA            NA          NA       NA        NA      NA
## MSQ4_AG_086          NA            NA          NA       NA        NA      NA
## MSQ4_AG_087          NA            NA          NA       NA        NA      NA
## MSQ4_AG_088          NA            NA          NA       NA        NA      NA
## MSQ4_AG_089          NA            NA          NA       NA        NA      NA
## MSQ4_AG_091          NA            NA          NA       NA        NA      NA
## MSQ4_AG_092          NA            NA          NA       NA        NA      NA
## MSQ4_AG_093          NA            NA          NA       NA        NA      NA
## MSQ4_AG_095          NA            NA          NA       NA        NA      NA
## MSQ4_AG_096          NA            NA          NA       NA        NA      NA
## MSQ4_AG_101          NA            NA          NA       NA        NA      NA
## MSQ4_AG_103          NA            NA          NA       NA        NA      NA
## MSQ4_AG_104          NA            NA          NA       NA        NA      NA
## MSQ4_AG_105          NA            NA          NA       NA        NA      NA
## MSQ4_AG_106          NA            NA          NA       NA        NA      NA
## MSQ4_AG_107          NA            NA          NA       NA        NA      NA
## MSQ4_AG_108          NA            NA          NA       NA        NA      NA
## MSQ4_AG_109          NA            NA          NA       NA        NA      NA
## MSQ4_AG_111          NA            NA          NA       NA        NA      NA
## MSQ4_AG_112          NA            NA          NA       NA        NA      NA
## MSQ4_AG_113          NA            NA          NA       NA        NA      NA
## MSQ4_AG_117          NA            NA          NA       NA        NA      NA
## MSQ4_AG_119          NA            NA          NA       NA        NA      NA
## MSQ4_AG_121          NA            NA          NA       NA        NA      NA
## MSQ4_AG_125          NA            NA          NA       NA        NA      NA
## MSQ4_AG_127          NA            NA          NA       NA        NA      NA
## MSQ4_AG_128          NA            NA          NA       NA        NA      NA
## MSQ4_AG_129          NA            NA          NA       NA        NA      NA
## MSQ4_AG_130          NA            NA          NA       NA        NA      NA
## MSQ4_AG_131          NA            NA          NA       NA        NA      NA
## MSQ4_AG_132          NA            NA          NA       NA        NA      NA
## MSQ4_AG_133          NA            NA          NA       NA        NA      NA
## MSQ4_AG_135          NA            NA          NA       NA        NA      NA
## MSQ4_AG_136          NA            NA          NA       NA        NA      NA
## MSQ4_AG_141          NA            NA          NA       NA        NA      NA
## MSQ4_AG_143          NA            NA          NA       NA        NA      NA
## MSQ4_AG_144          NA            NA          NA       NA        NA      NA
## MSQ4_AG_145          NA            NA          NA       NA        NA      NA
## MSQ4_AG_148          NA            NA          NA       NA        NA      NA
## MSQ4_AG_149          NA            NA          NA       NA        NA      NA
## MSQ4_AG_151          NA            NA          NA       NA        NA      NA
## MSQ4_AG_152          NA            NA          NA       NA        NA      NA
## MSQ4_AG_153          NA            NA          NA       NA        NA      NA
## MSQ4_AG_154          NA            NA          NA       NA        NA      NA
## MSQ4_AG_155          NA            NA          NA       NA        NA      NA
## MSQ4_AG_156          NA            NA          NA       NA        NA      NA
## MSQ4_AG_157          NA            NA          NA       NA        NA      NA
## MSQ4_AG_159          NA            NA          NA       NA        NA      NA
## MSQ4_AG_160          NA            NA          NA       NA        NA      NA
## MSQ4_AG_161          NA            NA          NA       NA        NA      NA
## MSQ4_AG_162          NA            NA          NA       NA        NA      NA
## MSQ4_AG_163          NA            NA          NA       NA        NA      NA
## MSQ4_AG_164          NA            NA          NA       NA        NA      NA
## MSQ4_AG_165          NA            NA          NA       NA        NA      NA
## MSQ4_AG_167          NA            NA          NA       NA        NA      NA
## MSQ4_AG_168          NA            NA          NA       NA        NA      NA
## MSQ4_AG_169          NA            NA          NA       NA        NA      NA
## MSQ4_AG_170          NA            NA          NA       NA        NA      NA
## MSQ4_AG_171          NA            NA          NA       NA        NA      NA
## MSQ4_AG_172          NA            NA          NA       NA        NA      NA
## MSQ4_AG_173          NA            NA          NA       NA        NA      NA
## MSQ4_AG_175          NA            NA          NA       NA        NA      NA
## MSQ4_AG_176          NA            NA          NA       NA        NA      NA
## MSQ4_AG_177          NA            NA          NA       NA        NA      NA
## MSQ4_AG_178          NA            NA          NA       NA        NA      NA
## MSQ4_AG_179          NA            NA          NA       NA        NA      NA
## MSQ4_AG_180          NA            NA          NA       NA        NA      NA
## MSQ4_AG_181          NA            NA          NA       NA        NA      NA
## MSQ4_AG_183          NA            NA          NA       NA        NA      NA
## MSQ4_AG_184          NA            NA          NA       NA        NA      NA
## MSQ4_AG_185          NA            NA          NA       NA        NA      NA
## MSQ4_AG_186          NA            NA          NA       NA        NA      NA
## MSQ4_AG_187          NA            NA          NA       NA        NA      NA
## MSQ4_AG_188          NA            NA          NA       NA        NA      NA
## MSQ4_AG_189          NA            NA          NA       NA        NA      NA
## MSQ4_AG_191          NA            NA          NA       NA        NA      NA
## MSQ4_AG_192          NA            NA          NA       NA        NA      NA
## MSQ4_AG_193          NA            NA          NA       NA        NA      NA
## MSQ4_AG_194          NA            NA          NA       NA        NA      NA
## MSQ4_AG_195          NA            NA          NA       NA        NA      NA
## MSQ4_AG_196          NA            NA          NA       NA        NA      NA
## MSQ4_AG_197          NA            NA          NA       NA        NA      NA
## MSQ4_AG_198          NA            NA          NA       NA        NA      NA
## MSQ4_AG_201          NA            NA          NA       NA        NA      NA
## MSQ4_AG_205          NA            NA          NA       NA        NA      NA
## MSQ4_AG_207          NA            NA          NA       NA        NA      NA
## MSQ4_AG_208          NA            NA          NA       NA        NA      NA
## MSQ4_AG_209          NA            NA          NA       NA        NA      NA
## MSQ4_AG_210          NA            NA          NA       NA        NA      NA
## MSQ4_AG_211          NA            NA          NA       NA        NA      NA
## MSQ4_AG_212          NA            NA          NA       NA        NA      NA
## MSQ4_AG_213          NA            NA          NA       NA        NA      NA
## MSQ4_AG_215          NA            NA          NA       NA        NA      NA
## MSQ4_AG_216          NA            NA          NA       NA        NA      NA
## MSQ4_AG_217          NA            NA          NA       NA        NA      NA
## MSQ4_AG_218          NA            NA          NA       NA        NA      NA
## MSQ4_AG_219          NA            NA          NA       NA        NA      NA
## MSQ4_AG_220          NA            NA          NA       NA        NA      NA
## MSQ4_AG_221          NA            NA          NA       NA        NA      NA
## MSQ4_AG_223          NA            NA          NA       NA        NA      NA
## MSQ4_AG_224          NA            NA          NA       NA        NA      NA
## MSQ4_AG_225          NA            NA          NA       NA        NA      NA
## MSQ4_AG_226          NA            NA          NA       NA        NA      NA
## MSQ4_AG_227          NA            NA          NA       NA        NA      NA
## MSQ4_AG_228          NA            NA          NA       NA        NA      NA
## MSQ4_AG_229          NA            NA          NA       NA        NA      NA
## MSQ4_AG_230          NA            NA          NA       NA        NA      NA
## MSQ4_AG_231          NA            NA          NA       NA        NA      NA
## MSQ4_AG_232          NA            NA          NA       NA        NA      NA
## MSQ4_AG_233          NA            NA          NA       NA        NA      NA
## MSQ4_AG_234          NA            NA          NA       NA        NA      NA
## MSQ4_AG_235          NA            NA          NA       NA        NA      NA
## MSQ4_AG_236          NA            NA          NA       NA        NA      NA
## MSQ4_AG_237          NA            NA          NA       NA        NA      NA
## MSQ4_AG_238          NA            NA          NA       NA        NA      NA
## MSQ4_AG_239          NA            NA          NA       NA        NA      NA
## MSQ4_AG_240          NA            NA          NA       NA        NA      NA
## MSQ4_AG_241          NA            NA          NA       NA        NA      NA
## MSQ4_AG_242          NA            NA          NA       NA        NA      NA
## MSQ4_AG_243          NA            NA          NA       NA        NA      NA
## MSQ4_AG_244          NA            NA          NA       NA        NA      NA
## MSQ4_AG_245          NA            NA          NA       NA        NA      NA
## MSQ4_AG_247          NA            NA          NA       NA        NA      NA
## MSQ4_AG_248          NA            NA          NA       NA        NA      NA
## MSQ4_AG_249          NA            NA          NA       NA        NA      NA
## MSQ6_AG_250          NA            NA          NA       NA        NA      NA
## MSQ6_AG_251          NA            NA          NA       NA        NA      NA
## MSQ6_AG_252          NA            NA          NA       NA        NA      NA
## MSQ6_AG_253          NA            NA          NA       NA        NA      NA
## MSQ6_AG_254          NA            NA          NA       NA        NA      NA
## MSQ6_AG_255          NA            NA          NA       NA        NA      NA
## MSQ6_AG_256          NA            NA          NA       NA        NA      NA
## MSQ6_AG_257          NA            NA          NA       NA        NA      NA
## MSQ6_AG_258          NA            NA          NA       NA        NA      NA
## MSQ6_AG_259          NA            NA          NA       NA        NA      NA
## MSQ6_AG_260          NA            NA          NA       NA        NA      NA
## MSQ6_AG_261          NA            NA          NA       NA        NA      NA
## MSQ6_AG_262          NA            NA          NA       NA        NA      NA
## MSQ6_AG_263          NA            NA          NA       NA        NA      NA
## MSQ6_AG_264          NA            NA          NA       NA        NA      NA
## MSQ6_AG_265          NA            NA          NA       NA        NA      NA
## MSQ6_AG_266          NA            NA          NA       NA        NA      NA
## MSQ6_AG_267          NA            NA          NA       NA        NA      NA
## MSQ6_AG_268          NA            NA          NA       NA        NA      NA
## MSQ6_AG_269          NA            NA          NA       NA        NA      NA
## MSQ6_AG_270          NA            NA          NA       NA        NA      NA
## MSQ6_AG_271          NA            NA          NA       NA        NA      NA
## MSQ6_AG_272          NA            NA          NA       NA        NA      NA
## MSQ6_AG_273          NA            NA          NA       NA        NA      NA
## MSQ6_AG_274          NA            NA          NA       NA        NA      NA
## MSQ6_AG_275          NA            NA          NA       NA        NA      NA
## MSQ6_AG_276          NA            NA          NA       NA        NA      NA
## MSQ6_AG_277          NA            NA          NA       NA        NA      NA
## MSQ6_AG_278          NA            NA          NA       NA        NA      NA
## MSQ6_AG_279          NA            NA          NA       NA        NA      NA
## MSQ6_AG_281          NA            NA          NA       NA        NA      NA
## MSQ6_AG_282          NA            NA          NA       NA        NA      NA
## MSQ6_AG_283          NA            NA          NA       NA        NA      NA
## MSQ6_AG_284          NA            NA          NA       NA        NA      NA
## MSQ6_AG_285          NA            NA          NA       NA        NA      NA
## MSQ6_AG_286          NA            NA          NA       NA        NA      NA
## MSQ6_AG_287          NA            NA          NA       NA        NA      NA
## MSQ6_AG_288          NA            NA          NA       NA        NA      NA
## MSQ6_AG_290          NA            NA          NA       NA        NA      NA
## MSQ6_AG_291          NA            NA          NA       NA        NA      NA
## MSQ6_AG_292          NA            NA          NA       NA        NA      NA
## MSQ6_AG_293          NA            NA          NA       NA        NA      NA
## MSQ6_AG_294          NA            NA          NA       NA        NA      NA
## MSQ6_AG_295          NA            NA          NA       NA        NA      NA
## MSQ6_AG_296          NA            NA          NA       NA        NA      NA
## MSQ6_AG_298          NA            NA          NA       NA        NA      NA
## MSQ6_AG_299          NA            NA          NA       NA        NA      NA
## MSQ6_AG_300          NA            NA          NA       NA        NA      NA
## MSQ6_AG_301          NA            NA          NA       NA        NA      NA
## MSQ6_AG_302          NA            NA          NA       NA        NA      NA
## MSQ6_AG_303          NA            NA          NA       NA        NA      NA
## MSQ6_AG_304          NA            NA          NA       NA        NA      NA
## MSQ6_AG_305          NA            NA          NA       NA        NA      NA
## MSQ6_AG_306          NA            NA          NA       NA        NA      NA
## MSQ6_AG_307          NA            NA          NA       NA        NA      NA
## MSQ6_AG_308          NA            NA          NA       NA        NA      NA
## MSQ6_AG_309          NA            NA          NA       NA        NA      NA
## MSQ6_AG_310          NA            NA          NA       NA        NA      NA
## MSQ6_AG_311          NA            NA          NA       NA        NA      NA
## MSQ6_AG_312          NA            NA          NA       NA        NA      NA
## MSQ6_AG_313          NA            NA          NA       NA        NA      NA
## MSQ6_AG_314          NA            NA          NA       NA        NA      NA
## MSQ6_AG_315          NA            NA          NA       NA        NA      NA
## MSQ6_AG_316          NA            NA          NA       NA        NA      NA
## MSQ6_AG_317          NA            NA          NA       NA        NA      NA
## MSQ6_AG_318          NA            NA          NA       NA        NA      NA
## MSQ6_AG_319          NA            NA          NA       NA        NA      NA
## MSQ6_AG_320          NA            NA          NA       NA        NA      NA
## MSQ6_AG_321          NA            NA          NA       NA        NA      NA
## MSQ6_AG_322          NA            NA          NA       NA        NA      NA
## MSQ6_AG_323          NA            NA          NA       NA        NA      NA
## MSQ6_AG_324          NA            NA          NA       NA        NA      NA
## MSQ6_AG_325          NA            NA          NA       NA        NA      NA
## MSQ6_AG_326          NA            NA          NA       NA        NA      NA
## MSQ6_AG_327          NA            NA          NA       NA        NA      NA
## MSQ6_AG_328          NA            NA          NA       NA        NA      NA
## MSQ6_AG_329          NA            NA          NA       NA        NA      NA
## MSQ6_AG_330          NA            NA          NA       NA        NA      NA
## MSQ6_AG_331          NA            NA          NA       NA        NA      NA
## MSQ6_AG_332          NA            NA          NA       NA        NA      NA
## MSQ6_AG_333          NA            NA          NA       NA        NA      NA
## MSQ6_AG_334          NA            NA          NA       NA        NA      NA
## MSQ6_AG_335          NA            NA          NA       NA        NA      NA
## MSQ6_AG_336          NA            NA          NA       NA        NA      NA
## MSQ6_AG_337          NA            NA          NA       NA        NA      NA
## MSQ6_AG_339          NA            NA          NA       NA        NA      NA
## MSQ6_AG_340          NA            NA          NA       NA        NA      NA
## MSQ6_AG_341          NA            NA          NA       NA        NA      NA
## MSQ6_AG_342          NA            NA          NA       NA        NA      NA
## MSQ6_AG_343          NA            NA          NA       NA        NA      NA
## MSQ6_AG_344          NA            NA          NA       NA        NA      NA
## MSQ6_AG_346          NA            NA          NA       NA        NA      NA
## MSQ6_AG_347          NA            NA          NA       NA        NA      NA
## MSQ6_AG_348          NA            NA          NA       NA        NA      NA
## MSQ6_AG_349          NA            NA          NA       NA        NA      NA
## MSQ6_AG_350          NA            NA          NA       NA        NA      NA
## MSQ6_AG_351          NA            NA          NA       NA        NA      NA
## MSQ6_AG_352          NA            NA          NA       NA        NA      NA
## MSQ6_AG_353          NA            NA          NA       NA        NA      NA
## MSQ6_AG_354          NA            NA          NA       NA        NA      NA
## MSQ6_AG_355          NA            NA          NA       NA        NA      NA
## MSQ6_AG_356          NA            NA          NA       NA        NA      NA
## MSQ6_AG_357          NA            NA          NA       NA        NA      NA
## MSQ6_AG_358          NA            NA          NA       NA        NA      NA
## MSQ6_AG_359          NA            NA          NA       NA        NA      NA
## MSQ6_AG_360          NA            NA          NA       NA        NA      NA
## MSQ6_AG_361          NA            NA          NA       NA        NA      NA
## MSQ6_AG_362          NA            NA          NA       NA        NA      NA
## MSQ6_AG_363          NA            NA          NA       NA        NA      NA
## MSQ6_AG_364          NA            NA          NA       NA        NA      NA
## MSQ6_AG_365          NA            NA          NA       NA        NA      NA
## MSQ6_AG_366          NA            NA          NA       NA        NA      NA
## MSQ6_AG_367          NA            NA          NA       NA        NA      NA
## MSQ6_AG_368          NA            NA          NA       NA        NA      NA
## MSQ6_AG_369          NA            NA          NA       NA        NA      NA
## MSQ6_AG_370          NA            NA          NA       NA        NA      NA
## MSQ6_AG_371          NA            NA          NA       NA        NA      NA
## MSQ6_AG_372          NA            NA          NA       NA        NA      NA
## MSQ6_AG_373          NA            NA          NA       NA        NA      NA
## MSQ6_AG_374          NA            NA          NA       NA        NA      NA
## MSQ6_AG_375          NA            NA          NA       NA        NA      NA
## MSQ6_AG_376          NA            NA          NA       NA        NA      NA
## MSQ6_AG_377          NA            NA          NA       NA        NA      NA
## MSQ6_AG_378          NA            NA          NA       NA        NA      NA
## MSQ6_AG_379          NA            NA          NA       NA        NA      NA
## MSQ6_AG_380          NA            NA          NA       NA        NA      NA
## MSQ6_AG_381          NA            NA          NA       NA        NA      NA
## MSQ6_AG_382          NA            NA          NA       NA        NA      NA
## MSQ6_AG_383          NA            NA          NA       NA        NA      NA
## MSQ6_AG_384          NA            NA          NA       NA        NA      NA
## MSQ6_AG_385          NA            NA          NA       NA        NA      NA
## MSQ6_AG_386          NA            NA          NA       NA        NA      NA
## MSQ6_AG_387          NA            NA          NA       NA        NA      NA
## MSQ6_AG_388      11.707         8.989      41.836    0.695        NA      NA
## MSQ6_AG_389      15.180         9.582      30.188    0.617        NA      NA
## MSQ6_AG_390      16.060         9.634      27.677    1.019        NA      NA
## MSQ6_AG_391      15.813        10.198      24.320    0.877        NA      NA
## MSQ6_AG_392       3.966         1.683      25.992    0.690        NA      NA
## MSQ6_AG_393      10.851         8.696      26.339    0.749        NA      NA
## MSQ6_AG_394      12.175        11.741      28.752    0.731        NA      NA
## MSQ6_AG_395       3.129         5.857      11.151    1.204        NA      NA
## MSQ6_AG_396      13.086         5.547      29.053    0.895        NA      NA
## MSQ6_AG_397      12.972         5.622      18.871    1.009        NA      NA
## MSQ6_AG_398      14.123         8.482      30.430    0.799        NA      NA
## MSQ6_AG_399      14.843         8.327      27.612    0.694        NA      NA
## MSQ6_AG_400      13.931         6.051      27.165    0.000        NA      NA
## MSQ6_AG_401      11.311         1.124      25.495    0.000        NA      NA
## MSQ6_AG_402      16.442         7.759      27.709    0.670        NA      NA
## MSQ6_AG_403      15.441         9.949      28.228    0.661        NA      NA
## MSQ6_AG_404      14.178         9.474      26.945    0.763        NA      NA
## MSQ6_AG_405      16.580         7.370      30.001    0.745        NA      NA
## MSQ6_AG_406      16.207         6.674      27.838    0.979        NA      NA
## MSQ6_AG_407      16.354         8.785      29.345    0.592        NA      NA
## MSQ6_AG_408      16.800         7.786      30.026    0.214        NA      NA
## MSQ6_AG_409      17.418         7.827      28.182    0.701        NA      NA
## MSQ6_AG_410      17.492        10.929      29.198    0.697        NA      NA
## MSQ6_AG_411      12.265         8.263      25.555    0.557        NA      NA
## MSQ6_AG_412      16.331         8.072      27.443    0.721        NA      NA
## MSQ6_AG_413      16.921         9.140      30.920    0.626        NA      NA
## MSQ6_AG_414      17.114        10.090      30.254    0.694        NA      NA
## MSQ6_AG_415      15.891         8.753      22.094    0.720        NA      NA
## MSQ6_AG_416      23.681        13.715      22.204    0.597        NA      NA
## MSQ6_AG_417      16.886         8.096      29.684    0.680        NA      NA
## MSQ6_AG_418      15.046         5.494      27.890    0.662        NA      NA
## MSQ6_AG_419      16.341         7.275      30.717    0.618        NA      NA
## MSQ6_AG_420      12.466         9.413      26.532    0.608        NA      NA
## MSQ6_AG_421      16.838         9.878      30.102    0.594        NA      NA
## MSQ6_AG_422      16.186         9.075      28.482    0.708        NA      NA
## MSQ6_AG_423      14.024         4.499      26.903    1.275        NA      NA
## MSQ6_AG_424      18.661         7.173      27.713    0.734        NA      NA
## MSQ6_AG_425      16.216         7.830      30.881    0.799        NA      NA
## MSQ6_AG_426      26.442         5.627      32.486    0.653        NA      NA
## MSQ6_AG_427      13.521         9.288      18.025    0.670        NA      NA
## MSQ6_AG_428      13.260        11.243      18.094    0.713        NA      NA
## MSQ6_AG_429       3.695         8.286      16.113    1.495        NA      NA
## MSQ6_AG_430       7.631         8.637      15.670    1.201        NA      NA
## MSQ6_AG_431       3.642         8.543      11.993    0.908        NA      NA
## MSQ6_AG_432      12.873         9.189      15.887    0.840        NA      NA
## MSQ6_AG_433      15.355         7.596      13.592    1.205        NA      NA
## MSQ6_AG_434      15.319         8.219      16.342    1.091        NA      NA
## MSQ6_AG_490          NA            NA          NA       NA        NA      NA
## MSQ6_AG_491          NA            NA          NA       NA        NA      NA
## MSQ6_AG_492          NA            NA          NA       NA        NA      NA
## MSQ6_AG_493          NA            NA          NA       NA        NA      NA
## MSQ6_AG_494          NA            NA          NA       NA        NA      NA
## MSQ6_AG_495          NA            NA          NA       NA        NA      NA
## MSQ6_AG_496          NA            NA          NA       NA        NA      NA
## MSQ6_AG_497          NA            NA          NA       NA        NA      NA
## MSQ6_AG_498          NA            NA          NA       NA        NA      NA
## MSQ6_AG_499          NA            NA          NA       NA        NA      NA
## MSQ6_AG_500          NA            NA          NA       NA        NA      NA
## MSQ6_AG_501          NA            NA          NA       NA        NA      NA
## MSQ6_AG_502          NA            NA          NA       NA        NA      NA
## MSQ6_AG_503          NA            NA          NA       NA        NA      NA
## MSQ6_AG_504          NA            NA          NA       NA        NA      NA
## MSQ6_AG_505          NA            NA          NA       NA        NA      NA
## MSQ6_AG_506          NA            NA          NA       NA        NA      NA
## MSQ6_AG_509          NA            NA          NA       NA        NA      NA
## MSQ6_AG_511          NA            NA          NA       NA        NA      NA
## MSQ6_AG_512          NA            NA          NA       NA        NA      NA
## MSQ6_AG_514          NA            NA          NA       NA        NA      NA
## MSQ6_AG_515          NA            NA          NA       NA        NA      NA
## MSQ6_AG_516          NA            NA          NA       NA        NA      NA
## MSQ6_AG_517          NA            NA          NA       NA        NA      NA
## MSQ6_AG_518          NA            NA          NA       NA        NA      NA
## MSQ6_AG_519          NA            NA          NA       NA        NA      NA
## MSQ6_AG_520          NA            NA          NA       NA        NA      NA
## MSQ6_AG_521          NA            NA          NA       NA        NA      NA
## MSQ6_AG_522          NA            NA          NA       NA        NA      NA
## MSQ6_AG_523          NA            NA          NA       NA        NA      NA
## MSQ6_AG_524          NA            NA          NA       NA        NA      NA
## MSQ6_AG_525          NA            NA          NA       NA        NA      NA
## MSQ6_AG_526          NA            NA          NA       NA        NA      NA
## MSQ6_AG_527          NA            NA          NA       NA        NA      NA
## MSQ6_AG_528          NA            NA          NA       NA        NA      NA
## MSQ6_AG_529          NA            NA          NA       NA        NA      NA
## MSQ6_AG_530          NA            NA          NA       NA        NA      NA
##             formate acetate propionate isobutyrate butyrate isovalerate
## MSQ4_AG_012      NA      NA         NA          NA       NA          NA
## MSQ4_AG_013      NA      NA         NA          NA       NA          NA
## MSQ4_AG_014      NA      NA         NA          NA       NA          NA
## MSQ4_AG_015      NA      NA         NA          NA       NA          NA
## MSQ4_AG_016      NA      NA         NA          NA       NA          NA
## MSQ4_AG_017      NA      NA         NA          NA       NA          NA
## MSQ4_AG_018      NA      NA         NA          NA       NA          NA
## MSQ4_AG_019      NA      NA         NA          NA       NA          NA
## MSQ4_AG_020      NA      NA         NA          NA       NA          NA
## MSQ4_AG_021      NA      NA         NA          NA       NA          NA
## MSQ4_AG_022      NA      NA         NA          NA       NA          NA
## MSQ4_AG_023      NA      NA         NA          NA       NA          NA
## MSQ4_AG_024      NA      NA         NA          NA       NA          NA
## MSQ4_AG_025      NA      NA         NA          NA       NA          NA
## MSQ4_AG_027      NA      NA         NA          NA       NA          NA
## MSQ4_AG_028      NA      NA         NA          NA       NA          NA
## MSQ4_AG_029      NA      NA         NA          NA       NA          NA
## MSQ4_AG_030      NA      NA         NA          NA       NA          NA
## MSQ4_AG_031      NA      NA         NA          NA       NA          NA
## MSQ4_AG_033      NA      NA         NA          NA       NA          NA
## MSQ4_AG_034      NA      NA         NA          NA       NA          NA
## MSQ4_AG_035      NA      NA         NA          NA       NA          NA
## MSQ4_AG_036      NA      NA         NA          NA       NA          NA
## MSQ4_AG_037      NA      NA         NA          NA       NA          NA
## MSQ4_AG_038      NA      NA         NA          NA       NA          NA
## MSQ4_AG_039      NA      NA         NA          NA       NA          NA
## MSQ4_AG_040      NA      NA         NA          NA       NA          NA
## MSQ4_AG_041      NA      NA         NA          NA       NA          NA
## MSQ4_AG_042      NA      NA         NA          NA       NA          NA
## MSQ4_AG_043      NA      NA         NA          NA       NA          NA
## MSQ4_AG_044      NA      NA         NA          NA       NA          NA
## MSQ4_AG_045      NA      NA         NA          NA       NA          NA
## MSQ4_AG_046      NA      NA         NA          NA       NA          NA
## MSQ4_AG_047      NA      NA         NA          NA       NA          NA
## MSQ4_AG_049      NA      NA         NA          NA       NA          NA
## MSQ4_AG_050      NA      NA         NA          NA       NA          NA
## MSQ4_AG_051      NA      NA         NA          NA       NA          NA
## MSQ4_AG_052      NA      NA         NA          NA       NA          NA
## MSQ4_AG_053      NA      NA         NA          NA       NA          NA
## MSQ4_AG_054      NA      NA         NA          NA       NA          NA
## MSQ4_AG_055      NA      NA         NA          NA       NA          NA
## MSQ4_AG_056      NA      NA         NA          NA       NA          NA
## MSQ4_AG_057      NA      NA         NA          NA       NA          NA
## MSQ4_AG_058      NA      NA         NA          NA       NA          NA
## MSQ4_AG_059      NA      NA         NA          NA       NA          NA
## MSQ4_AG_062      NA      NA         NA          NA       NA          NA
## MSQ4_AG_064      NA      NA         NA          NA       NA          NA
## MSQ4_AG_065      NA      NA         NA          NA       NA          NA
## MSQ4_AG_072      NA      NA         NA          NA       NA          NA
## MSQ4_AG_076      NA      NA         NA          NA       NA          NA
## MSQ4_AG_077      NA      NA         NA          NA       NA          NA
## MSQ4_AG_080      NA      NA         NA          NA       NA          NA
## MSQ4_AG_083      NA      NA         NA          NA       NA          NA
## MSQ4_AG_084      NA      NA         NA          NA       NA          NA
## MSQ4_AG_085      NA      NA         NA          NA       NA          NA
## MSQ4_AG_086      NA      NA         NA          NA       NA          NA
## MSQ4_AG_087      NA      NA         NA          NA       NA          NA
## MSQ4_AG_088      NA      NA         NA          NA       NA          NA
## MSQ4_AG_089      NA      NA         NA          NA       NA          NA
## MSQ4_AG_091      NA      NA         NA          NA       NA          NA
## MSQ4_AG_092      NA      NA         NA          NA       NA          NA
## MSQ4_AG_093      NA      NA         NA          NA       NA          NA
## MSQ4_AG_095      NA      NA         NA          NA       NA          NA
## MSQ4_AG_096      NA      NA         NA          NA       NA          NA
## MSQ4_AG_101      NA      NA         NA          NA       NA          NA
## MSQ4_AG_103      NA      NA         NA          NA       NA          NA
## MSQ4_AG_104      NA      NA         NA          NA       NA          NA
## MSQ4_AG_105      NA      NA         NA          NA       NA          NA
## MSQ4_AG_106      NA      NA         NA          NA       NA          NA
## MSQ4_AG_107      NA      NA         NA          NA       NA          NA
## MSQ4_AG_108      NA      NA         NA          NA       NA          NA
## MSQ4_AG_109      NA      NA         NA          NA       NA          NA
## MSQ4_AG_111      NA      NA         NA          NA       NA          NA
## MSQ4_AG_112      NA      NA         NA          NA       NA          NA
## MSQ4_AG_113      NA      NA         NA          NA       NA          NA
## MSQ4_AG_117      NA      NA         NA          NA       NA          NA
## MSQ4_AG_119      NA      NA         NA          NA       NA          NA
## MSQ4_AG_121      NA      NA         NA          NA       NA          NA
## MSQ4_AG_125      NA      NA         NA          NA       NA          NA
## MSQ4_AG_127      NA      NA         NA          NA       NA          NA
## MSQ4_AG_128      NA      NA         NA          NA       NA          NA
## MSQ4_AG_129      NA      NA         NA          NA       NA          NA
## MSQ4_AG_130      NA      NA         NA          NA       NA          NA
## MSQ4_AG_131      NA      NA         NA          NA       NA          NA
## MSQ4_AG_132      NA      NA         NA          NA       NA          NA
## MSQ4_AG_133      NA      NA         NA          NA       NA          NA
## MSQ4_AG_135      NA      NA         NA          NA       NA          NA
## MSQ4_AG_136      NA      NA         NA          NA       NA          NA
## MSQ4_AG_141      NA      NA         NA          NA       NA          NA
## MSQ4_AG_143      NA      NA         NA          NA       NA          NA
## MSQ4_AG_144      NA      NA         NA          NA       NA          NA
## MSQ4_AG_145      NA      NA         NA          NA       NA          NA
## MSQ4_AG_148      NA      NA         NA          NA       NA          NA
## MSQ4_AG_149      NA      NA         NA          NA       NA          NA
## MSQ4_AG_151      NA      NA         NA          NA       NA          NA
## MSQ4_AG_152      NA      NA         NA          NA       NA          NA
## MSQ4_AG_153      NA      NA         NA          NA       NA          NA
## MSQ4_AG_154      NA      NA         NA          NA       NA          NA
## MSQ4_AG_155      NA      NA         NA          NA       NA          NA
## MSQ4_AG_156      NA      NA         NA          NA       NA          NA
## MSQ4_AG_157      NA      NA         NA          NA       NA          NA
## MSQ4_AG_159      NA      NA         NA          NA       NA          NA
## MSQ4_AG_160      NA      NA         NA          NA       NA          NA
## MSQ4_AG_161      NA      NA         NA          NA       NA          NA
## MSQ4_AG_162      NA      NA         NA          NA       NA          NA
## MSQ4_AG_163      NA      NA         NA          NA       NA          NA
## MSQ4_AG_164      NA      NA         NA          NA       NA          NA
## MSQ4_AG_165      NA      NA         NA          NA       NA          NA
## MSQ4_AG_167      NA      NA         NA          NA       NA          NA
## MSQ4_AG_168      NA      NA         NA          NA       NA          NA
## MSQ4_AG_169      NA      NA         NA          NA       NA          NA
## MSQ4_AG_170      NA      NA         NA          NA       NA          NA
## MSQ4_AG_171      NA      NA         NA          NA       NA          NA
## MSQ4_AG_172      NA      NA         NA          NA       NA          NA
## MSQ4_AG_173      NA      NA         NA          NA       NA          NA
## MSQ4_AG_175      NA      NA         NA          NA       NA          NA
## MSQ4_AG_176      NA      NA         NA          NA       NA          NA
## MSQ4_AG_177      NA      NA         NA          NA       NA          NA
## MSQ4_AG_178      NA      NA         NA          NA       NA          NA
## MSQ4_AG_179      NA      NA         NA          NA       NA          NA
## MSQ4_AG_180      NA      NA         NA          NA       NA          NA
## MSQ4_AG_181      NA      NA         NA          NA       NA          NA
## MSQ4_AG_183      NA      NA         NA          NA       NA          NA
## MSQ4_AG_184      NA      NA         NA          NA       NA          NA
## MSQ4_AG_185      NA      NA         NA          NA       NA          NA
## MSQ4_AG_186      NA      NA         NA          NA       NA          NA
## MSQ4_AG_187      NA      NA         NA          NA       NA          NA
## MSQ4_AG_188      NA      NA         NA          NA       NA          NA
## MSQ4_AG_189      NA      NA         NA          NA       NA          NA
## MSQ4_AG_191      NA      NA         NA          NA       NA          NA
## MSQ4_AG_192      NA      NA         NA          NA       NA          NA
## MSQ4_AG_193      NA      NA         NA          NA       NA          NA
## MSQ4_AG_194      NA      NA         NA          NA       NA          NA
## MSQ4_AG_195      NA      NA         NA          NA       NA          NA
## MSQ4_AG_196      NA      NA         NA          NA       NA          NA
## MSQ4_AG_197      NA      NA         NA          NA       NA          NA
## MSQ4_AG_198      NA      NA         NA          NA       NA          NA
## MSQ4_AG_201      NA      NA         NA          NA       NA          NA
## MSQ4_AG_205      NA      NA         NA          NA       NA          NA
## MSQ4_AG_207      NA      NA         NA          NA       NA          NA
## MSQ4_AG_208      NA      NA         NA          NA       NA          NA
## MSQ4_AG_209      NA      NA         NA          NA       NA          NA
## MSQ4_AG_210      NA      NA         NA          NA       NA          NA
## MSQ4_AG_211      NA      NA         NA          NA       NA          NA
## MSQ4_AG_212      NA      NA         NA          NA       NA          NA
## MSQ4_AG_213      NA      NA         NA          NA       NA          NA
## MSQ4_AG_215      NA      NA         NA          NA       NA          NA
## MSQ4_AG_216      NA      NA         NA          NA       NA          NA
## MSQ4_AG_217      NA      NA         NA          NA       NA          NA
## MSQ4_AG_218      NA      NA         NA          NA       NA          NA
## MSQ4_AG_219      NA      NA         NA          NA       NA          NA
## MSQ4_AG_220      NA      NA         NA          NA       NA          NA
## MSQ4_AG_221      NA      NA         NA          NA       NA          NA
## MSQ4_AG_223      NA      NA         NA          NA       NA          NA
## MSQ4_AG_224      NA      NA         NA          NA       NA          NA
## MSQ4_AG_225      NA      NA         NA          NA       NA          NA
## MSQ4_AG_226      NA      NA         NA          NA       NA          NA
## MSQ4_AG_227      NA      NA         NA          NA       NA          NA
## MSQ4_AG_228      NA      NA         NA          NA       NA          NA
## MSQ4_AG_229      NA      NA         NA          NA       NA          NA
## MSQ4_AG_230      NA      NA         NA          NA       NA          NA
## MSQ4_AG_231      NA      NA         NA          NA       NA          NA
## MSQ4_AG_232      NA      NA         NA          NA       NA          NA
## MSQ4_AG_233      NA      NA         NA          NA       NA          NA
## MSQ4_AG_234      NA      NA         NA          NA       NA          NA
## MSQ4_AG_235      NA      NA         NA          NA       NA          NA
## MSQ4_AG_236      NA      NA         NA          NA       NA          NA
## MSQ4_AG_237      NA      NA         NA          NA       NA          NA
## MSQ4_AG_238      NA      NA         NA          NA       NA          NA
## MSQ4_AG_239      NA      NA         NA          NA       NA          NA
## MSQ4_AG_240      NA      NA         NA          NA       NA          NA
## MSQ4_AG_241      NA      NA         NA          NA       NA          NA
## MSQ4_AG_242      NA      NA         NA          NA       NA          NA
## MSQ4_AG_243      NA      NA         NA          NA       NA          NA
## MSQ4_AG_244      NA      NA         NA          NA       NA          NA
## MSQ4_AG_245      NA      NA         NA          NA       NA          NA
## MSQ4_AG_247      NA      NA         NA          NA       NA          NA
## MSQ4_AG_248      NA      NA         NA          NA       NA          NA
## MSQ4_AG_249      NA      NA         NA          NA       NA          NA
## MSQ6_AG_250      NA      NA         NA          NA       NA          NA
## MSQ6_AG_251      NA      NA         NA          NA       NA          NA
## MSQ6_AG_252      NA      NA         NA          NA       NA          NA
## MSQ6_AG_253      NA      NA         NA          NA       NA          NA
## MSQ6_AG_254      NA      NA         NA          NA       NA          NA
## MSQ6_AG_255      NA      NA         NA          NA       NA          NA
## MSQ6_AG_256      NA      NA         NA          NA       NA          NA
## MSQ6_AG_257      NA      NA         NA          NA       NA          NA
## MSQ6_AG_258      NA      NA         NA          NA       NA          NA
## MSQ6_AG_259      NA      NA         NA          NA       NA          NA
## MSQ6_AG_260      NA      NA         NA          NA       NA          NA
## MSQ6_AG_261      NA      NA         NA          NA       NA          NA
## MSQ6_AG_262      NA      NA         NA          NA       NA          NA
## MSQ6_AG_263      NA      NA         NA          NA       NA          NA
## MSQ6_AG_264      NA      NA         NA          NA       NA          NA
## MSQ6_AG_265      NA      NA         NA          NA       NA          NA
## MSQ6_AG_266      NA      NA         NA          NA       NA          NA
## MSQ6_AG_267      NA      NA         NA          NA       NA          NA
## MSQ6_AG_268      NA      NA         NA          NA       NA          NA
## MSQ6_AG_269      NA      NA         NA          NA       NA          NA
## MSQ6_AG_270      NA      NA         NA          NA       NA          NA
## MSQ6_AG_271      NA      NA         NA          NA       NA          NA
## MSQ6_AG_272      NA      NA         NA          NA       NA          NA
## MSQ6_AG_273      NA      NA         NA          NA       NA          NA
## MSQ6_AG_274      NA      NA         NA          NA       NA          NA
## MSQ6_AG_275      NA      NA         NA          NA       NA          NA
## MSQ6_AG_276      NA      NA         NA          NA       NA          NA
## MSQ6_AG_277      NA      NA         NA          NA       NA          NA
## MSQ6_AG_278      NA      NA         NA          NA       NA          NA
## MSQ6_AG_279      NA      NA         NA          NA       NA          NA
## MSQ6_AG_281      NA      NA         NA          NA       NA          NA
## MSQ6_AG_282      NA      NA         NA          NA       NA          NA
## MSQ6_AG_283      NA      NA         NA          NA       NA          NA
## MSQ6_AG_284      NA      NA         NA          NA       NA          NA
## MSQ6_AG_285      NA      NA         NA          NA       NA          NA
## MSQ6_AG_286      NA      NA         NA          NA       NA          NA
## MSQ6_AG_287      NA      NA         NA          NA       NA          NA
## MSQ6_AG_288      NA      NA         NA          NA       NA          NA
## MSQ6_AG_290      NA      NA         NA          NA       NA          NA
## MSQ6_AG_291      NA      NA         NA          NA       NA          NA
## MSQ6_AG_292      NA      NA         NA          NA       NA          NA
## MSQ6_AG_293      NA      NA         NA          NA       NA          NA
## MSQ6_AG_294      NA      NA         NA          NA       NA          NA
## MSQ6_AG_295      NA      NA         NA          NA       NA          NA
## MSQ6_AG_296      NA      NA         NA          NA       NA          NA
## MSQ6_AG_298      NA      NA         NA          NA       NA          NA
## MSQ6_AG_299      NA      NA         NA          NA       NA          NA
## MSQ6_AG_300      NA      NA         NA          NA       NA          NA
## MSQ6_AG_301      NA      NA         NA          NA       NA          NA
## MSQ6_AG_302      NA      NA         NA          NA       NA          NA
## MSQ6_AG_303      NA      NA         NA          NA       NA          NA
## MSQ6_AG_304      NA      NA         NA          NA       NA          NA
## MSQ6_AG_305      NA      NA         NA          NA       NA          NA
## MSQ6_AG_306      NA      NA         NA          NA       NA          NA
## MSQ6_AG_307      NA      NA         NA          NA       NA          NA
## MSQ6_AG_308      NA      NA         NA          NA       NA          NA
## MSQ6_AG_309      NA      NA         NA          NA       NA          NA
## MSQ6_AG_310      NA      NA         NA          NA       NA          NA
## MSQ6_AG_311      NA      NA         NA          NA       NA          NA
## MSQ6_AG_312      NA      NA         NA          NA       NA          NA
## MSQ6_AG_313      NA      NA         NA          NA       NA          NA
## MSQ6_AG_314      NA      NA         NA          NA       NA          NA
## MSQ6_AG_315      NA      NA         NA          NA       NA          NA
## MSQ6_AG_316      NA      NA         NA          NA       NA          NA
## MSQ6_AG_317      NA      NA         NA          NA       NA          NA
## MSQ6_AG_318      NA      NA         NA          NA       NA          NA
## MSQ6_AG_319      NA      NA         NA          NA       NA          NA
## MSQ6_AG_320      NA      NA         NA          NA       NA          NA
## MSQ6_AG_321      NA      NA         NA          NA       NA          NA
## MSQ6_AG_322      NA      NA         NA          NA       NA          NA
## MSQ6_AG_323      NA      NA         NA          NA       NA          NA
## MSQ6_AG_324      NA      NA         NA          NA       NA          NA
## MSQ6_AG_325      NA      NA         NA          NA       NA          NA
## MSQ6_AG_326      NA      NA         NA          NA       NA          NA
## MSQ6_AG_327      NA      NA         NA          NA       NA          NA
## MSQ6_AG_328      NA      NA         NA          NA       NA          NA
## MSQ6_AG_329      NA      NA         NA          NA       NA          NA
## MSQ6_AG_330      NA      NA         NA          NA       NA          NA
## MSQ6_AG_331      NA      NA         NA          NA       NA          NA
## MSQ6_AG_332      NA      NA         NA          NA       NA          NA
## MSQ6_AG_333      NA      NA         NA          NA       NA          NA
## MSQ6_AG_334      NA      NA         NA          NA       NA          NA
## MSQ6_AG_335      NA      NA         NA          NA       NA          NA
## MSQ6_AG_336      NA      NA         NA          NA       NA          NA
## MSQ6_AG_337      NA      NA         NA          NA       NA          NA
## MSQ6_AG_339      NA      NA         NA          NA       NA          NA
## MSQ6_AG_340      NA      NA         NA          NA       NA          NA
## MSQ6_AG_341      NA      NA         NA          NA       NA          NA
## MSQ6_AG_342      NA      NA         NA          NA       NA          NA
## MSQ6_AG_343      NA      NA         NA          NA       NA          NA
## MSQ6_AG_344      NA      NA         NA          NA       NA          NA
## MSQ6_AG_346      NA      NA         NA          NA       NA          NA
## MSQ6_AG_347      NA      NA         NA          NA       NA          NA
## MSQ6_AG_348      NA      NA         NA          NA       NA          NA
## MSQ6_AG_349      NA      NA         NA          NA       NA          NA
## MSQ6_AG_350      NA      NA         NA          NA       NA          NA
## MSQ6_AG_351      NA      NA         NA          NA       NA          NA
## MSQ6_AG_352      NA      NA         NA          NA       NA          NA
## MSQ6_AG_353      NA      NA         NA          NA       NA          NA
## MSQ6_AG_354      NA      NA         NA          NA       NA          NA
## MSQ6_AG_355      NA      NA         NA          NA       NA          NA
## MSQ6_AG_356      NA      NA         NA          NA       NA          NA
## MSQ6_AG_357      NA      NA         NA          NA       NA          NA
## MSQ6_AG_358      NA      NA         NA          NA       NA          NA
## MSQ6_AG_359      NA      NA         NA          NA       NA          NA
## MSQ6_AG_360      NA      NA         NA          NA       NA          NA
## MSQ6_AG_361      NA      NA         NA          NA       NA          NA
## MSQ6_AG_362      NA      NA         NA          NA       NA          NA
## MSQ6_AG_363      NA      NA         NA          NA       NA          NA
## MSQ6_AG_364      NA      NA         NA          NA       NA          NA
## MSQ6_AG_365      NA      NA         NA          NA       NA          NA
## MSQ6_AG_366      NA      NA         NA          NA       NA          NA
## MSQ6_AG_367      NA      NA         NA          NA       NA          NA
## MSQ6_AG_368      NA      NA         NA          NA       NA          NA
## MSQ6_AG_369      NA      NA         NA          NA       NA          NA
## MSQ6_AG_370      NA      NA         NA          NA       NA          NA
## MSQ6_AG_371      NA      NA         NA          NA       NA          NA
## MSQ6_AG_372      NA      NA         NA          NA       NA          NA
## MSQ6_AG_373      NA      NA         NA          NA       NA          NA
## MSQ6_AG_374      NA      NA         NA          NA       NA          NA
## MSQ6_AG_375      NA      NA         NA          NA       NA          NA
## MSQ6_AG_376      NA      NA         NA          NA       NA          NA
## MSQ6_AG_377      NA      NA         NA          NA       NA          NA
## MSQ6_AG_378      NA      NA         NA          NA       NA          NA
## MSQ6_AG_379      NA      NA         NA          NA       NA          NA
## MSQ6_AG_380      NA      NA         NA          NA       NA          NA
## MSQ6_AG_381      NA      NA         NA          NA       NA          NA
## MSQ6_AG_382      NA      NA         NA          NA       NA          NA
## MSQ6_AG_383      NA      NA         NA          NA       NA          NA
## MSQ6_AG_384      NA      NA         NA          NA       NA          NA
## MSQ6_AG_385      NA      NA         NA          NA       NA          NA
## MSQ6_AG_386      NA      NA         NA          NA       NA          NA
## MSQ6_AG_387      NA      NA         NA          NA       NA          NA
## MSQ6_AG_388      NA      NA         NA          NA       NA          NA
## MSQ6_AG_389      NA      NA         NA          NA       NA          NA
## MSQ6_AG_390      NA      NA         NA          NA       NA          NA
## MSQ6_AG_391      NA      NA         NA          NA       NA          NA
## MSQ6_AG_392      NA      NA         NA          NA       NA          NA
## MSQ6_AG_393      NA      NA         NA          NA       NA          NA
## MSQ6_AG_394      NA      NA         NA          NA       NA          NA
## MSQ6_AG_395      NA      NA         NA          NA       NA          NA
## MSQ6_AG_396      NA      NA         NA          NA       NA          NA
## MSQ6_AG_397      NA      NA         NA          NA       NA          NA
## MSQ6_AG_398      NA      NA         NA          NA       NA          NA
## MSQ6_AG_399      NA      NA         NA          NA       NA          NA
## MSQ6_AG_400      NA      NA         NA          NA       NA          NA
## MSQ6_AG_401      NA      NA         NA          NA       NA          NA
## MSQ6_AG_402      NA      NA         NA          NA       NA          NA
## MSQ6_AG_403      NA      NA         NA          NA       NA          NA
## MSQ6_AG_404      NA      NA         NA          NA       NA          NA
## MSQ6_AG_405      NA      NA         NA          NA       NA          NA
## MSQ6_AG_406      NA      NA         NA          NA       NA          NA
## MSQ6_AG_407      NA      NA         NA          NA       NA          NA
## MSQ6_AG_408      NA      NA         NA          NA       NA          NA
## MSQ6_AG_409      NA      NA         NA          NA       NA          NA
## MSQ6_AG_410      NA      NA         NA          NA       NA          NA
## MSQ6_AG_411      NA      NA         NA          NA       NA          NA
## MSQ6_AG_412      NA      NA         NA          NA       NA          NA
## MSQ6_AG_413      NA      NA         NA          NA       NA          NA
## MSQ6_AG_414      NA      NA         NA          NA       NA          NA
## MSQ6_AG_415      NA      NA         NA          NA       NA          NA
## MSQ6_AG_416      NA      NA         NA          NA       NA          NA
## MSQ6_AG_417      NA      NA         NA          NA       NA          NA
## MSQ6_AG_418      NA      NA         NA          NA       NA          NA
## MSQ6_AG_419      NA      NA         NA          NA       NA          NA
## MSQ6_AG_420      NA      NA         NA          NA       NA          NA
## MSQ6_AG_421      NA      NA         NA          NA       NA          NA
## MSQ6_AG_422      NA      NA         NA          NA       NA          NA
## MSQ6_AG_423      NA      NA         NA          NA       NA          NA
## MSQ6_AG_424      NA      NA         NA          NA       NA          NA
## MSQ6_AG_425      NA      NA         NA          NA       NA          NA
## MSQ6_AG_426      NA      NA         NA          NA       NA          NA
## MSQ6_AG_427      NA      NA         NA          NA       NA          NA
## MSQ6_AG_428      NA      NA         NA          NA       NA          NA
## MSQ6_AG_429      NA      NA         NA          NA       NA          NA
## MSQ6_AG_430      NA      NA         NA          NA       NA          NA
## MSQ6_AG_431      NA      NA         NA          NA       NA          NA
## MSQ6_AG_432      NA      NA         NA          NA       NA          NA
## MSQ6_AG_433      NA      NA         NA          NA       NA          NA
## MSQ6_AG_434      NA      NA         NA          NA       NA          NA
## MSQ6_AG_490      NA      NA         NA          NA       NA          NA
## MSQ6_AG_491      NA      NA         NA          NA       NA          NA
## MSQ6_AG_492      NA      NA         NA          NA       NA          NA
## MSQ6_AG_493      NA      NA         NA          NA       NA          NA
## MSQ6_AG_494      NA      NA         NA          NA       NA          NA
## MSQ6_AG_495      NA      NA         NA          NA       NA          NA
## MSQ6_AG_496      NA      NA         NA          NA       NA          NA
## MSQ6_AG_497      NA      NA         NA          NA       NA          NA
## MSQ6_AG_498      NA      NA         NA          NA       NA          NA
## MSQ6_AG_499      NA      NA         NA          NA       NA          NA
## MSQ6_AG_500      NA      NA         NA          NA       NA          NA
## MSQ6_AG_501      NA      NA         NA          NA       NA          NA
## MSQ6_AG_502      NA      NA         NA          NA       NA          NA
## MSQ6_AG_503      NA      NA         NA          NA       NA          NA
## MSQ6_AG_504      NA      NA         NA          NA       NA          NA
## MSQ6_AG_505      NA      NA         NA          NA       NA          NA
## MSQ6_AG_506      NA      NA         NA          NA       NA          NA
## MSQ6_AG_509      NA      NA         NA          NA       NA          NA
## MSQ6_AG_511      NA      NA         NA          NA       NA          NA
## MSQ6_AG_512      NA      NA         NA          NA       NA          NA
## MSQ6_AG_514      NA      NA         NA          NA       NA          NA
## MSQ6_AG_515      NA      NA         NA          NA       NA          NA
## MSQ6_AG_516      NA      NA         NA          NA       NA          NA
## MSQ6_AG_517      NA      NA         NA          NA       NA          NA
## MSQ6_AG_518      NA      NA         NA          NA       NA          NA
## MSQ6_AG_519      NA      NA         NA          NA       NA          NA
## MSQ6_AG_520      NA      NA         NA          NA       NA          NA
## MSQ6_AG_521      NA      NA         NA          NA       NA          NA
## MSQ6_AG_522      NA      NA         NA          NA       NA          NA
## MSQ6_AG_523      NA      NA         NA          NA       NA          NA
## MSQ6_AG_524      NA      NA         NA          NA       NA          NA
## MSQ6_AG_525      NA      NA         NA          NA       NA          NA
## MSQ6_AG_526      NA      NA         NA          NA       NA          NA
## MSQ6_AG_527      NA      NA         NA          NA       NA          NA
## MSQ6_AG_528      NA      NA         NA          NA       NA          NA
## MSQ6_AG_529      NA      NA         NA          NA       NA          NA
## MSQ6_AG_530      NA      NA         NA          NA       NA          NA
##             valerate BCFA total is.neg period day_fermentation_d IR_TR
## MSQ4_AG_012       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_013       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_014       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_015       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_016       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_017       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_018       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_019       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_020       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_021       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_022       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_023       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_024       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_025       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_027       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_028       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_029       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_030       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_031       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_033       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_034       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_035       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_036       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_037       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_038       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_039       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_040       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_041       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_042       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_043       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_044       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_045       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_046       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_047       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_049       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_050       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_051       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_052       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_053       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_054       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_055       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_056       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_057       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_058       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_059       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_062       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_064       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_065       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_072       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_076       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_077       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_080       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_083       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_084       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_085       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_086       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_087       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_088       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_089       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_091       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_092       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_093       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_095       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_096       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_101       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_103       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_104       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_105       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_106       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_107       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_108       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_109       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_111       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_112       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_113       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_117       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_119       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_121       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_125       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_127       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_128       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_129       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_130       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_131       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_132       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_133       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_135       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_136       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_141       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_143       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_144       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_145       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_148       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_149       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_151       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_152       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_153       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_154       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_155       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_156       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_157       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_159       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_160       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_161       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_162       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_163       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_164       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_165       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_167       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_168       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_169       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_170       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_171       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_172       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_173       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_175       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_176       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_177       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_178       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_179       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_180       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_181       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_183       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_184       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_185       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_186       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_187       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_188       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_189       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_191       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_192       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_193       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_194       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_195       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_196       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_197       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_198       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_201       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_205       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_207       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_208       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_209       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_210       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_211       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_212       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_213       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_215       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_216       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_217       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_218       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_219       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_220       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_221       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_223       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_224       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_225       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_226       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_227       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_228       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_229       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_230       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_231       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_232       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_233       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_234       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_235       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_236       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_237       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_238       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_239       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_240       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_241       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_242       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_243       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_244       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_245       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_247       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_248       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ4_AG_249       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_250       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_251       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_252       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_253       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_254       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_255       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_256       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_257       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_258       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_259       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_260       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_261       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_262       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_263       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_264       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_265       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_266       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_267       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_268       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_269       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_270       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_271       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_272       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_273       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_274       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_275       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_276       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_277       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_278       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_279       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_281       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_282       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_283       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_284       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_285       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_286       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_287       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_288       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_290       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_291       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_292       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_293       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_294       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_295       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_296       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_298       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_299       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_300       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_301       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_302       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_303       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_304       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_305       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_306       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_307       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_308       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_309       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_310       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_311       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_312       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_313       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_314       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_315       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_316       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_317       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_318       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_319       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_320       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_321       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_322       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_323       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_324       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_325       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_326       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_327       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_328       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_329       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_330       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_331       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_332       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_333       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_334       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_335       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_336       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_337       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_339       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_340       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_341       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_342       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_343       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_344       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_346       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_347       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_348       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_349       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_350       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_351       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_352       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_353       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_354       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_355       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_356       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_357       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_358       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_359       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_360       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_361       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_362       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_363       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_364       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_365       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_366       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_367       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_368       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_369       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_370       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_371       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_372       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_373       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_374       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_375       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_376       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_377       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_378       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_379       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_380       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_381       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_382       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_383       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_384       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_385       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_386       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_387       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_388       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_389       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_390       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_391       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_392       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_393       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_394       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_395       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_396       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_397       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_398       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_399       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_400       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_401       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_402       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_403       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_404       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_405       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_406       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_407       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_408       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_409       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_410       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_411       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_412       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_413       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_414       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_415       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_416       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_417       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_418       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_419       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_420       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_421       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_422       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_423       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_424       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_425       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_426       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_427       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_428       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_429       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_430       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_431       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_432       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_433       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_434       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_490       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_491       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_492       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_493       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_494       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_495       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_496       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_497       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_498       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_499       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_500       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_501       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_502       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_503       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_504       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_505       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_506       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_509       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_511       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_512       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_514       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_515       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_516       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_517       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_518       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_519       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_520       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_521       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_522       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_523       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_524       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_525       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_526       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_527       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_528       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_529       NA   NA    NA  FALSE    pNA                 NA  <NA>
## MSQ6_AG_530       NA   NA    NA  FALSE    pNA                 NA  <NA>
##             Stab_Treat experiment_period period_treat period_reactor
## MSQ4_AG_012       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_013       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_014       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_015       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_016       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_017       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_018       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_019       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_020       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_021       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_022       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_023       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_024       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_025       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_027       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_028       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_029       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_030       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_031       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_033       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_034       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_035       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_036       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_037       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_038       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_039       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_040       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_041       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_042       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_043       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_044       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_045       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_046       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_047       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_049       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_050       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_051       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_052       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_053       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_054       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_055       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_056       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_057       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_058       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_059       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_062       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_064       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_065       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_072       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_076       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_077       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_080       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_083       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_084       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_085       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_086       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_087       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_088       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_089       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_091       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_092       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_093       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_095       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_096       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_101       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_103       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_104       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_105       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_106       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_107       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_108       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_109       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_111       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_112       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_113       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_117       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_119       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_121       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_125       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_127       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_128       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_129       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_130       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_131       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_132       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_133       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_135       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_136       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_141       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_143       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_144       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_145       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_148       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_149       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_151       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_152       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_153       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_154       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_155       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_156       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_157       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_159       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_160       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_161       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_162       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_163       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_164       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_165       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_167       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_168       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_169       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_170       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_171       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_172       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_173       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_175       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_176       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_177       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_178       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_179       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_180       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_181       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_183       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_184       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_185       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_186       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_187       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_188       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_189       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_191       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_192       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_193       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_194       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_195       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_196       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_197       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_198       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_201       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_205       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_207       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_208       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_209       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_210       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_211       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_212       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_213       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_215       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_216       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_217       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_218       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_219       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_220       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_221       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_223       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_224       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_225       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_226       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_227       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_228       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_229       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_230       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_231       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_232       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_233       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_234       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_235       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_236       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_237       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_238       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_239       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_240       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_241       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_242       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_243       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_244       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_245       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_247       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_248       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ4_AG_249       <NA>         mIMT_1_NA         <NA>         pNA_NA
## MSQ6_AG_250       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_251       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_252       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_253       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_254       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_255       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_256       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_257       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_258       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_259       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_260       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_261       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_262       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_263       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_264       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_265       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_266       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_267       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_268       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_269       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_270       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_271       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_272       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_273       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_274       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_275       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_276       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_277       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_278       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_279       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_281       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_282       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_283       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_284       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_285       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_286       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_287       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_288       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_290       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_291       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_292       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_293       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_294       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_295       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_296       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_298       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_299       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_300       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_301       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_302       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_303       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_304       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_305       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_306       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_307       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_308       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_309       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_310       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_311       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_312       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_313       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_314       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_315       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_316       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_317       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_318       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_319       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_320       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_321       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_322       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_323       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_324       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_325       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_326       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_327       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_328       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_329       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_330       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_331       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_332       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_333       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_334       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_335       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_336       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_337       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_339       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_340       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_341       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_342       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_343       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_344       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_346       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_347       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_348       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_349       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_350       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_351       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_352       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_353       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_354       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_355       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_356       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_357       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_358       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_359       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_360       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_361       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_362       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_363       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_364       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_365       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_366       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_367       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_368       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_369       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_370       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_371       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_372       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_373       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_374       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_375       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_376       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_377       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_378       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_379       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_380       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_381       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_382       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_383       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_384       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_385       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_386       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_387       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_388       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_389       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_390       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_391       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_392       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_393       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_394       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_395       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_396       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_397       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_398       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_399       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_400       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_401       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_402       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_403       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_404       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_405       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_406       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_407       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_408       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_409       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_410       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_411       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_412       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_413       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_414       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_415       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_416       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_417       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_418       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_419       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_420       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_421       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_422       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_423       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_424       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_425       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_426       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_427       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_428       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_429       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_430       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_431       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_432       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_433       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_434       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_490       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_491       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_492       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_493       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_494       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_495       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_496       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_497       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_498       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_499       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_500       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_501       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_502       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_503       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_504       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_505       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_506       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_509       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_511       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_512       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_514       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_515       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_516       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_517       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_518       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_519       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_520       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_521       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_522       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_523       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_524       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_525       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_526       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_527       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_528       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_529       <NA>         mIMT_2_NA         <NA>         pNA_NA
## MSQ6_AG_530       <NA>         mIMT_2_NA         <NA>         pNA_NA
##             treatment_grouped treatment_invivo day_num
## MSQ4_AG_012               H2O             none      -1
## MSQ4_AG_013               H2O             none      -1
## MSQ4_AG_014               H2O             none      -1
## MSQ4_AG_015               H2O             none      -1
## MSQ4_AG_016               H2O             none      -1
## MSQ4_AG_017               H2O             none      -1
## MSQ4_AG_018               H2O             none      -1
## MSQ4_AG_019               H2O             none      -1
## MSQ4_AG_020               H2O             none      -1
## MSQ4_AG_021               H2O             none      -1
## MSQ4_AG_022               DSS             none      -1
## MSQ4_AG_023               DSS             none      -1
## MSQ4_AG_024               DSS             none      -1
## MSQ4_AG_025               DSS             none      -1
## MSQ4_AG_027               DSS             none      -1
## MSQ4_AG_028               DSS             none      -1
## MSQ4_AG_029               DSS             none      -1
## MSQ4_AG_030               DSS             none      -1
## MSQ4_AG_031               DSS             none      -1
## MSQ4_AG_033          Eubiotic             none      -1
## MSQ4_AG_034          Eubiotic             none      -1
## MSQ4_AG_035          Eubiotic             none      -1
## MSQ4_AG_036          Eubiotic             none      -1
## MSQ4_AG_037         Dysbiotic             none      -1
## MSQ4_AG_038         Dysbiotic             none      -1
## MSQ4_AG_039         Dysbiotic             none      -1
## MSQ4_AG_040         Dysbiotic             none      -1
## MSQ4_AG_041         Dysbiotic             none      -1
## MSQ4_AG_042       Normobiotic             none      -1
## MSQ4_AG_043       Normobiotic             none      -1
## MSQ4_AG_044       Normobiotic             none      -1
## MSQ4_AG_045       Normobiotic             none      -1
## MSQ4_AG_046       Normobiotic             none      -1
## MSQ4_AG_047               H2O             none       6
## MSQ4_AG_049               H2O             none       6
## MSQ4_AG_050               H2O             none       6
## MSQ4_AG_051               H2O             none       6
## MSQ4_AG_052               H2O             none       6
## MSQ4_AG_053               H2O             none       6
## MSQ4_AG_054               H2O             none       6
## MSQ4_AG_055               H2O             none       6
## MSQ4_AG_056               H2O             none       6
## MSQ4_AG_057               DSS              DSS       6
## MSQ4_AG_058               DSS              DSS       6
## MSQ4_AG_059               DSS              DSS       6
## MSQ4_AG_062               DSS              DSS       6
## MSQ4_AG_064               DSS              DSS       6
## MSQ4_AG_065               DSS              DSS       6
## MSQ4_AG_072         Dysbiotic        Dysbiotic       6
## MSQ4_AG_076         Dysbiotic        Dysbiotic       6
## MSQ4_AG_077       Normobiotic      Normobiotic       6
## MSQ4_AG_080       Normobiotic      Normobiotic       6
## MSQ4_AG_083               H2O             none      12
## MSQ4_AG_084               H2O             none      12
## MSQ4_AG_085               H2O             none      12
## MSQ4_AG_086               H2O             none      12
## MSQ4_AG_087               H2O             none      12
## MSQ4_AG_088               H2O             none      12
## MSQ4_AG_089               H2O             none      12
## MSQ4_AG_091               H2O             none      12
## MSQ4_AG_092               DSS              DSS      12
## MSQ4_AG_093               DSS              DSS      12
## MSQ4_AG_095               DSS              DSS      12
## MSQ4_AG_096               DSS              DSS      12
## MSQ4_AG_101          Eubiotic         Eubiotic      12
## MSQ4_AG_103          Eubiotic         Eubiotic      12
## MSQ4_AG_104          Eubiotic         Eubiotic      12
## MSQ4_AG_105          Eubiotic         Eubiotic      12
## MSQ4_AG_106         Dysbiotic        Dysbiotic      12
## MSQ4_AG_107         Dysbiotic        Dysbiotic      12
## MSQ4_AG_108         Dysbiotic        Dysbiotic      12
## MSQ4_AG_109         Dysbiotic        Dysbiotic      12
## MSQ4_AG_111       Normobiotic      Normobiotic      12
## MSQ4_AG_112       Normobiotic      Normobiotic      12
## MSQ4_AG_113       Normobiotic      Normobiotic      12
## MSQ4_AG_117               H2O             none      14
## MSQ4_AG_119               H2O             none      14
## MSQ4_AG_121               H2O             none      14
## MSQ4_AG_125               H2O             none      14
## MSQ4_AG_127               DSS              DSS      14
## MSQ4_AG_128               DSS              DSS      14
## MSQ4_AG_129               DSS              DSS      14
## MSQ4_AG_130               DSS              DSS      14
## MSQ4_AG_131               DSS              DSS      14
## MSQ4_AG_132               DSS              DSS      14
## MSQ4_AG_133               DSS              DSS      14
## MSQ4_AG_135          Eubiotic         Eubiotic      14
## MSQ4_AG_136          Eubiotic         Eubiotic      14
## MSQ4_AG_141         Dysbiotic        Dysbiotic      14
## MSQ4_AG_143         Dysbiotic        Dysbiotic      14
## MSQ4_AG_144         Dysbiotic        Dysbiotic      14
## MSQ4_AG_145       Normobiotic      Normobiotic      14
## MSQ4_AG_148       Normobiotic      Normobiotic      14
## MSQ4_AG_149       Normobiotic      Normobiotic      14
## MSQ4_AG_151               H2O             none      14
## MSQ4_AG_152               H2O             none      14
## MSQ4_AG_153               H2O             none      14
## MSQ4_AG_154               H2O             none      14
## MSQ4_AG_155               H2O             none      14
## MSQ4_AG_156               H2O             none      14
## MSQ4_AG_157               H2O             none      14
## MSQ4_AG_159               H2O             none      14
## MSQ4_AG_160               DSS              DSS      14
## MSQ4_AG_161               DSS              DSS      14
## MSQ4_AG_162               DSS              DSS      14
## MSQ4_AG_163               DSS              DSS      14
## MSQ4_AG_164               DSS              DSS      14
## MSQ4_AG_165               DSS              DSS      14
## MSQ4_AG_167               DSS              DSS      14
## MSQ4_AG_168               DSS              DSS      14
## MSQ4_AG_169          Eubiotic         Eubiotic      14
## MSQ4_AG_170          Eubiotic         Eubiotic      14
## MSQ4_AG_171          Eubiotic         Eubiotic      14
## MSQ4_AG_172          Eubiotic         Eubiotic      14
## MSQ4_AG_173          Eubiotic         Eubiotic      14
## MSQ4_AG_175         Dysbiotic        Dysbiotic      14
## MSQ4_AG_176         Dysbiotic        Dysbiotic      14
## MSQ4_AG_177         Dysbiotic        Dysbiotic      14
## MSQ4_AG_178         Dysbiotic        Dysbiotic      14
## MSQ4_AG_179         Dysbiotic        Dysbiotic      14
## MSQ4_AG_180       Normobiotic      Normobiotic      14
## MSQ4_AG_181       Normobiotic      Normobiotic      14
## MSQ4_AG_183       Normobiotic      Normobiotic      14
## MSQ4_AG_184               H2O             none      14
## MSQ4_AG_185               H2O             none      14
## MSQ4_AG_186               H2O             none      14
## MSQ4_AG_187               H2O             none      14
## MSQ4_AG_188               H2O             none      14
## MSQ4_AG_189               H2O             none      14
## MSQ4_AG_191               H2O             none      14
## MSQ4_AG_192               H2O             none      14
## MSQ4_AG_193               H2O             none      14
## MSQ4_AG_194               DSS              DSS      14
## MSQ4_AG_195               DSS              DSS      14
## MSQ4_AG_196               DSS              DSS      14
## MSQ4_AG_197               DSS              DSS      14
## MSQ4_AG_198               DSS              DSS      14
## MSQ4_AG_201               DSS              DSS      14
## MSQ4_AG_205          Eubiotic         Eubiotic      14
## MSQ4_AG_207         Dysbiotic        Dysbiotic      14
## MSQ4_AG_208         Dysbiotic        Dysbiotic      14
## MSQ4_AG_209         Dysbiotic        Dysbiotic      14
## MSQ4_AG_210       Normobiotic      Normobiotic      14
## MSQ4_AG_211       Normobiotic      Normobiotic      14
## MSQ4_AG_212       Normobiotic      Normobiotic      14
## MSQ4_AG_213       Normobiotic      Normobiotic      14
## MSQ4_AG_215               H2O             none      12
## MSQ4_AG_216               H2O             none      12
## MSQ4_AG_217               H2O             none      12
## MSQ4_AG_218               H2O             none      12
## MSQ4_AG_219               H2O             none      12
## MSQ4_AG_220               H2O             none      12
## MSQ4_AG_221               H2O             none      12
## MSQ4_AG_223               DSS              DSS      12
## MSQ4_AG_224               DSS              DSS      12
## MSQ4_AG_225               DSS              DSS      12
## MSQ4_AG_226               DSS              DSS      12
## MSQ4_AG_227               DSS              DSS      12
## MSQ4_AG_228               DSS              DSS      12
## MSQ4_AG_229          Eubiotic         Eubiotic      12
## MSQ4_AG_230          Eubiotic         Eubiotic      12
## MSQ4_AG_231          Eubiotic         Eubiotic      12
## MSQ4_AG_232          Eubiotic         Eubiotic      12
## MSQ4_AG_233         Dysbiotic        Dysbiotic      12
## MSQ4_AG_234         Dysbiotic        Dysbiotic      12
## MSQ4_AG_235         Dysbiotic        Dysbiotic      12
## MSQ4_AG_236         Dysbiotic        Dysbiotic      12
## MSQ4_AG_237       Normobiotic      Normobiotic      12
## MSQ4_AG_238       Normobiotic      Normobiotic      12
## MSQ4_AG_239       Normobiotic      Normobiotic      12
## MSQ4_AG_240       Normobiotic      Normobiotic      12
## MSQ4_AG_241               H2O             none      -1
## MSQ4_AG_242               H2O             none      -1
## MSQ4_AG_243               H2O             none      -1
## MSQ4_AG_244               DSS             none      -1
## MSQ4_AG_245               DSS             none      -1
## MSQ4_AG_247         Dysbiotic             none      -1
## MSQ4_AG_248       Normobiotic             none      -1
## MSQ4_AG_249       Normobiotic             none      -1
## MSQ6_AG_250          Eubiotic             none      -1
## MSQ6_AG_251          Eubiotic             none      -1
## MSQ6_AG_252          Eubiotic             none      -1
## MSQ6_AG_253          Eubiotic             none      -1
## MSQ6_AG_254          Eubiotic             none      -1
## MSQ6_AG_255          Eubiotic             none      -1
## MSQ6_AG_256          Eubiotic             none      -1
## MSQ6_AG_257          Eubiotic             none      -1
## MSQ6_AG_258         Dysbiotic             none      -1
## MSQ6_AG_259         Dysbiotic             none      -1
## MSQ6_AG_260         Dysbiotic             none      -1
## MSQ6_AG_261         Dysbiotic             none      -1
## MSQ6_AG_262         Dysbiotic             none      -1
## MSQ6_AG_263         Dysbiotic             none      -1
## MSQ6_AG_264         Dysbiotic             none      -1
## MSQ6_AG_265         Dysbiotic             none      -1
## MSQ6_AG_266       Normobiotic             none      -1
## MSQ6_AG_267       Normobiotic             none      -1
## MSQ6_AG_268       Normobiotic             none      -1
## MSQ6_AG_269       Normobiotic             none      -1
## MSQ6_AG_270       Normobiotic             none      -1
## MSQ6_AG_271       Normobiotic             none      -1
## MSQ6_AG_272       Normobiotic             none      -1
## MSQ6_AG_273       Normobiotic             none      -1
## MSQ6_AG_274            Predni             none      -1
## MSQ6_AG_275            Predni             none      -1
## MSQ6_AG_276            Predni             none      -1
## MSQ6_AG_277            Predni             none      -1
## MSQ6_AG_278            Predni             none      -1
## MSQ6_AG_279            Predni             none      -1
## MSQ6_AG_281            Predni             none      -1
## MSQ6_AG_282               DSS             none      -1
## MSQ6_AG_283               DSS             none      -1
## MSQ6_AG_284               DSS             none      -1
## MSQ6_AG_285               DSS             none      -1
## MSQ6_AG_286               DSS             none      -1
## MSQ6_AG_287               DSS             none      -1
## MSQ6_AG_288               DSS             none      -1
## MSQ6_AG_290               H2O             none      -1
## MSQ6_AG_291               H2O             none      -1
## MSQ6_AG_292               H2O             none      -1
## MSQ6_AG_293               H2O             none      -1
## MSQ6_AG_294               H2O             none      -1
## MSQ6_AG_295               H2O             none      -1
## MSQ6_AG_296               H2O             none      -1
## MSQ6_AG_298          Eubiotic         Eubiotic      12
## MSQ6_AG_299          Eubiotic         Eubiotic      12
## MSQ6_AG_300          Eubiotic         Eubiotic      12
## MSQ6_AG_301          Eubiotic         Eubiotic      12
## MSQ6_AG_302          Eubiotic         Eubiotic      12
## MSQ6_AG_303          Eubiotic         Eubiotic      12
## MSQ6_AG_304          Eubiotic         Eubiotic      12
## MSQ6_AG_305         Dysbiotic        Dysbiotic      12
## MSQ6_AG_306         Dysbiotic        Dysbiotic      12
## MSQ6_AG_307         Dysbiotic        Dysbiotic      12
## MSQ6_AG_308         Dysbiotic        Dysbiotic      12
## MSQ6_AG_309         Dysbiotic        Dysbiotic      12
## MSQ6_AG_310         Dysbiotic        Dysbiotic      12
## MSQ6_AG_311         Dysbiotic        Dysbiotic      12
## MSQ6_AG_312         Dysbiotic        Dysbiotic      12
## MSQ6_AG_313       Normobiotic      Normobiotic      12
## MSQ6_AG_314       Normobiotic      Normobiotic      12
## MSQ6_AG_315       Normobiotic      Normobiotic      12
## MSQ6_AG_316       Normobiotic      Normobiotic      12
## MSQ6_AG_317       Normobiotic      Normobiotic      12
## MSQ6_AG_318       Normobiotic      Normobiotic      12
## MSQ6_AG_319       Normobiotic      Normobiotic      12
## MSQ6_AG_320       Normobiotic      Normobiotic      12
## MSQ6_AG_321            Predni           Predni      12
## MSQ6_AG_322            Predni           Predni      12
## MSQ6_AG_323            Predni           Predni      12
## MSQ6_AG_324            Predni           Predni      12
## MSQ6_AG_325            Predni           Predni      12
## MSQ6_AG_326            Predni           Predni      12
## MSQ6_AG_327            Predni           Predni      12
## MSQ6_AG_328            Predni           Predni      12
## MSQ6_AG_329               DSS              DSS      12
## MSQ6_AG_330               DSS              DSS      12
## MSQ6_AG_331               DSS              DSS      12
## MSQ6_AG_332               DSS              DSS      12
## MSQ6_AG_333               DSS              DSS      12
## MSQ6_AG_334               DSS              DSS      12
## MSQ6_AG_335               DSS              DSS      12
## MSQ6_AG_336               H2O             none      12
## MSQ6_AG_337               H2O             none      12
## MSQ6_AG_339               H2O             none      12
## MSQ6_AG_340               H2O             none      12
## MSQ6_AG_341               H2O             none      12
## MSQ6_AG_342               H2O             none      12
## MSQ6_AG_343               H2O             none      12
## MSQ6_AG_344          Eubiotic         Eubiotic      14
## MSQ6_AG_346          Eubiotic         Eubiotic      14
## MSQ6_AG_347          Eubiotic         Eubiotic      14
## MSQ6_AG_348          Eubiotic         Eubiotic      14
## MSQ6_AG_349          Eubiotic         Eubiotic      14
## MSQ6_AG_350          Eubiotic         Eubiotic      14
## MSQ6_AG_351         Dysbiotic        Dysbiotic      14
## MSQ6_AG_352         Dysbiotic        Dysbiotic      14
## MSQ6_AG_353         Dysbiotic        Dysbiotic      14
## MSQ6_AG_354         Dysbiotic        Dysbiotic      14
## MSQ6_AG_355         Dysbiotic        Dysbiotic      14
## MSQ6_AG_356         Dysbiotic        Dysbiotic      14
## MSQ6_AG_357         Dysbiotic        Dysbiotic      14
## MSQ6_AG_358         Dysbiotic        Dysbiotic      14
## MSQ6_AG_359       Normobiotic      Normobiotic      14
## MSQ6_AG_360       Normobiotic      Normobiotic      14
## MSQ6_AG_361       Normobiotic      Normobiotic      14
## MSQ6_AG_362       Normobiotic      Normobiotic      14
## MSQ6_AG_363       Normobiotic      Normobiotic      14
## MSQ6_AG_364       Normobiotic      Normobiotic      14
## MSQ6_AG_365            Predni           Predni      14
## MSQ6_AG_366            Predni           Predni      14
## MSQ6_AG_367            Predni           Predni      14
## MSQ6_AG_368            Predni           Predni      14
## MSQ6_AG_369            Predni           Predni      14
## MSQ6_AG_370            Predni           Predni      14
## MSQ6_AG_371            Predni           Predni      14
## MSQ6_AG_372               DSS              DSS      14
## MSQ6_AG_373               DSS              DSS      14
## MSQ6_AG_374               DSS              DSS      14
## MSQ6_AG_375               DSS              DSS      14
## MSQ6_AG_376               DSS              DSS      14
## MSQ6_AG_377               DSS              DSS      14
## MSQ6_AG_378               DSS              DSS      14
## MSQ6_AG_379               DSS              DSS      14
## MSQ6_AG_380               H2O             none      14
## MSQ6_AG_381               H2O             none      14
## MSQ6_AG_382               H2O             none      14
## MSQ6_AG_383               H2O             none      14
## MSQ6_AG_384               H2O             none      14
## MSQ6_AG_385               H2O             none      14
## MSQ6_AG_386               H2O             none      14
## MSQ6_AG_387               H2O             none      14
## MSQ6_AG_388          Eubiotic         Eubiotic      14
## MSQ6_AG_389          Eubiotic         Eubiotic      14
## MSQ6_AG_390          Eubiotic         Eubiotic      14
## MSQ6_AG_391          Eubiotic         Eubiotic      14
## MSQ6_AG_392          Eubiotic         Eubiotic      14
## MSQ6_AG_393          Eubiotic         Eubiotic      14
## MSQ6_AG_394          Eubiotic         Eubiotic      14
## MSQ6_AG_395         Dysbiotic        Dysbiotic      14
## MSQ6_AG_396         Dysbiotic        Dysbiotic      14
## MSQ6_AG_397         Dysbiotic        Dysbiotic      14
## MSQ6_AG_398         Dysbiotic        Dysbiotic      14
## MSQ6_AG_399         Dysbiotic        Dysbiotic      14
## MSQ6_AG_400         Dysbiotic        Dysbiotic      14
## MSQ6_AG_401         Dysbiotic        Dysbiotic      14
## MSQ6_AG_402         Dysbiotic        Dysbiotic      14
## MSQ6_AG_403       Normobiotic      Normobiotic      14
## MSQ6_AG_404       Normobiotic      Normobiotic      14
## MSQ6_AG_405       Normobiotic      Normobiotic      14
## MSQ6_AG_406       Normobiotic      Normobiotic      14
## MSQ6_AG_407       Normobiotic      Normobiotic      14
## MSQ6_AG_408       Normobiotic      Normobiotic      14
## MSQ6_AG_409       Normobiotic      Normobiotic      14
## MSQ6_AG_410       Normobiotic      Normobiotic      14
## MSQ6_AG_411            Predni           Predni      14
## MSQ6_AG_412            Predni           Predni      14
## MSQ6_AG_413            Predni           Predni      14
## MSQ6_AG_414            Predni           Predni      14
## MSQ6_AG_415            Predni           Predni      14
## MSQ6_AG_416            Predni           Predni      14
## MSQ6_AG_417            Predni           Predni      14
## MSQ6_AG_418            Predni           Predni      14
## MSQ6_AG_419               DSS              DSS      14
## MSQ6_AG_420               DSS              DSS      14
## MSQ6_AG_421               DSS              DSS      14
## MSQ6_AG_422               DSS              DSS      14
## MSQ6_AG_423               DSS              DSS      14
## MSQ6_AG_424               DSS              DSS      14
## MSQ6_AG_425               DSS              DSS      14
## MSQ6_AG_426               DSS              DSS      14
## MSQ6_AG_427               H2O             none      14
## MSQ6_AG_428               H2O             none      14
## MSQ6_AG_429               H2O             none      14
## MSQ6_AG_430               H2O             none      14
## MSQ6_AG_431               H2O             none      14
## MSQ6_AG_432               H2O             none      14
## MSQ6_AG_433               H2O             none      14
## MSQ6_AG_434               H2O             none      14
## MSQ6_AG_490          Eubiotic         Eubiotic       6
## MSQ6_AG_491          Eubiotic         Eubiotic       6
## MSQ6_AG_492          Eubiotic         Eubiotic       6
## MSQ6_AG_493          Eubiotic         Eubiotic       6
## MSQ6_AG_494          Eubiotic         Eubiotic       6
## MSQ6_AG_495         Dysbiotic        Dysbiotic       6
## MSQ6_AG_496         Dysbiotic        Dysbiotic       6
## MSQ6_AG_497         Dysbiotic        Dysbiotic       6
## MSQ6_AG_498         Dysbiotic        Dysbiotic       6
## MSQ6_AG_499         Dysbiotic        Dysbiotic       6
## MSQ6_AG_500         Dysbiotic        Dysbiotic       6
## MSQ6_AG_501         Dysbiotic        Dysbiotic       6
## MSQ6_AG_502         Dysbiotic        Dysbiotic       6
## MSQ6_AG_503       Normobiotic      Normobiotic       6
## MSQ6_AG_504       Normobiotic      Normobiotic       6
## MSQ6_AG_505       Normobiotic      Normobiotic       6
## MSQ6_AG_506       Normobiotic      Normobiotic       6
## MSQ6_AG_509       Normobiotic      Normobiotic       6
## MSQ6_AG_511            Predni           Predni       6
## MSQ6_AG_512            Predni           Predni       6
## MSQ6_AG_514            Predni           Predni       6
## MSQ6_AG_515            Predni           Predni       6
## MSQ6_AG_516            Predni           Predni       6
## MSQ6_AG_517               DSS              DSS       6
## MSQ6_AG_518               DSS              DSS       6
## MSQ6_AG_519               DSS              DSS       6
## MSQ6_AG_520               DSS              DSS       6
## MSQ6_AG_521               DSS              DSS       6
## MSQ6_AG_522               DSS              DSS       6
## MSQ6_AG_523               H2O             none       6
## MSQ6_AG_524               H2O             none       6
## MSQ6_AG_525               H2O             none       6
## MSQ6_AG_526               H2O             none       6
## MSQ6_AG_527               H2O             none       6
## MSQ6_AG_528               H2O             none       6
## MSQ6_AG_529               H2O             none       6
## MSQ6_AG_530               H2O             none       6
```


```r
physeq %>% 
  subset_samples(sample_type == "Feces" & 
                   experiment == "mIMT_2" &
                   treatment_grouped %in% c("H2O", 
                                            "Eubiotic",
                                            "Dysbiotic",
                                            "DSS")) -> ps_invivo

sample_names(ps_invivo) <- sample_data(ps_invivo)$Sample


ps_invivo %>% 
  sample_data() %>% 
  data.frame()
```

```
##       sample_file  input filtered filtered_pc denoisedF denoisedR denoisedF_pc
## S_219 MSQ6_AG_250  44852    44233        0.99     42985     43681         0.97
## S_220 MSQ6_AG_251  72545    72048        0.99     70726     70818         0.98
## S_221 MSQ6_AG_252  76560    76051        0.99     74773     75012         0.98
## S_222 MSQ6_AG_253  62479    62192        1.00     61528     61201         0.99
## S_223 MSQ6_AG_254  78944    78670        1.00     77717     77643         0.99
## S_224 MSQ6_AG_255  72930    72622        1.00     71294     71600         0.98
## S_225 MSQ6_AG_256 194545   193572        0.99    191373    191090         0.99
## S_226 MSQ6_AG_257  62034    61405        0.99     60035     60599         0.98
## S_227 MSQ6_AG_258  77255    76602        0.99     75092     75664         0.98
## S_228 MSQ6_AG_259  58784    58444        0.99     57469     57504         0.98
## S_229 MSQ6_AG_260  65923    65701        1.00     64821     65102         0.99
## S_230 MSQ6_AG_261  51065    50870        1.00     50385     50108         0.99
## S_231 MSQ6_AG_262 123400   123134        1.00    122255    122251         0.99
## S_232 MSQ6_AG_263  60612    60294        0.99     59542     59544         0.99
## S_233 MSQ6_AG_264  66974    66745        1.00     65951     65893         0.99
## S_234 MSQ6_AG_265  59225    58828        0.99     57814     57973         0.98
## S_251 MSQ6_AG_282  46134    45771        0.99     44977     45491         0.98
## S_252 MSQ6_AG_283  43200    42935        0.99     41926     42309         0.98
## S_253 MSQ6_AG_284  61771    61392        0.99     60360     60789         0.98
## S_254 MSQ6_AG_285  52064    51842        1.00     51293     51367         0.99
## S_255 MSQ6_AG_286  60985    60720        1.00     60048     60142         0.99
## S_256 MSQ6_AG_287  58725    58433        1.00     57619     57826         0.99
## S_257 MSQ6_AG_288  53218    52990        1.00     52344     52408         0.99
## S_259 MSQ6_AG_290  35191    34713        0.99     33684     34288         0.97
## S_260 MSQ6_AG_291  59786    59402        0.99     58312     58638         0.98
## S_261 MSQ6_AG_292  71389    70960        0.99     69815     70122         0.98
## S_262 MSQ6_AG_293  48071    47899        1.00     47393     47356         0.99
## S_263 MSQ6_AG_294  60459    60267        1.00     59414     59486         0.99
## S_264 MSQ6_AG_295  56699    56368        0.99     55608     55822         0.99
## S_265 MSQ6_AG_296  52415    52164        1.00     51488     51508         0.99
## S_267 MSQ6_AG_298  40629    40017        0.98     37962     39383         0.95
## S_268 MSQ6_AG_299  67522    67104        0.99     65052     65909         0.97
## S_269 MSQ6_AG_300  62759    62323        0.99     60618     61516         0.97
## S_270 MSQ6_AG_301  51923    51673        1.00     50793     51026         0.98
## S_271 MSQ6_AG_302  67574    67308        1.00     66518     66567         0.99
## S_272 MSQ6_AG_303  69843    69434        0.99     68180     68564         0.98
## S_273 MSQ6_AG_304  53162    52908        1.00     52171     52277         0.99
## S_274 MSQ6_AG_305  89399    88860        0.99     86777     87744         0.98
## S_275 MSQ6_AG_306  47024    46664        0.99     46031     46150         0.99
## S_276 MSQ6_AG_307  52772    52491        0.99     51730     51714         0.99
## S_277 MSQ6_AG_308  57920    57707        1.00     56978     56801         0.99
## S_278 MSQ6_AG_309  45769    45629        1.00     45263     44921         0.99
## S_279 MSQ6_AG_310  69581    69393        1.00     68667     68504         0.99
## S_280 MSQ6_AG_311  51654    51454        1.00     50897     50887         0.99
## S_281 MSQ6_AG_312  47303    47147        1.00     46707     46496         0.99
## S_298 MSQ6_AG_329  45576    45320        0.99     44563     44714         0.98
## S_299 MSQ6_AG_330  37551    37094        0.99     36164     36676         0.97
## S_300 MSQ6_AG_331  49047    48750        0.99     47852     48161         0.98
## S_301 MSQ6_AG_332  55629    55307        0.99     54626     54653         0.99
## S_302 MSQ6_AG_333  39112    38954        1.00     38658     38534         0.99
## S_303 MSQ6_AG_334  56158    55984        1.00     54978     55559         0.98
## S_304 MSQ6_AG_335  56174    55990        1.00     54558     55582         0.97
## S_305 MSQ6_AG_336  49080    48860        1.00     48261     48353         0.99
## S_306 MSQ6_AG_337  51145    50825        0.99     50082     50314         0.99
## S_308 MSQ6_AG_339  61065    60799        1.00     59746     59875         0.98
## S_309 MSQ6_AG_340  49604    49210        0.99     48411     48539         0.98
## S_310 MSQ6_AG_341  39676    39568        1.00     39079     39070         0.99
## S_311 MSQ6_AG_342  67509    67341        1.00     66594     66566         0.99
## S_312 MSQ6_AG_343  54562    54327        1.00     53556     53619         0.99
## S_313 MSQ6_AG_344  51990    51826        1.00     51343     51218         0.99
## S_315 MSQ6_AG_346  56389    56103        0.99     55311     55282         0.99
## S_316 MSQ6_AG_347  52602    52284        0.99     51592     51797         0.99
## S_317 MSQ6_AG_348  55353    55045        0.99     54304     54490         0.99
## S_318 MSQ6_AG_349  47228    47005        1.00     46369     46415         0.99
## S_319 MSQ6_AG_350  51051    50825        1.00     50302     50300         0.99
## S_320 MSQ6_AG_351  50339    50129        1.00     49527     49544         0.99
## S_321 MSQ6_AG_352  53955    53733        1.00     53035     53076         0.99
## S_322 MSQ6_AG_353  50561    50357        1.00     49892     49982         0.99
## S_323 MSQ6_AG_354  55667    55424        1.00     54749     54623         0.99
## S_324 MSQ6_AG_355  53070    52783        0.99     52078     52396         0.99
## S_325 MSQ6_AG_356  43764    43485        0.99     42777     43048         0.98
## S_326 MSQ6_AG_357  55187    55049        1.00     54210     54461         0.98
## S_327 MSQ6_AG_358  44766    44572        1.00     44063     44115         0.99
## S_341 MSQ6_AG_372  51618    51384        1.00     50750     50912         0.99
## S_342 MSQ6_AG_373  38695    38576        1.00     38190     38209         0.99
## S_343 MSQ6_AG_374  31639    31476        0.99     31038     31086         0.99
## S_344 MSQ6_AG_375  57314    57167        1.00     56524     56611         0.99
## S_345 MSQ6_AG_376  26796    26719        1.00     26323     26354         0.99
## S_346 MSQ6_AG_377  42867    42762        1.00     42236     42312         0.99
## S_347 MSQ6_AG_378  49591    49374        1.00     48708     48943         0.99
## S_348 MSQ6_AG_379  49706    49595        1.00     49100     49181         0.99
## S_349 MSQ6_AG_380  82068    81706        1.00     80665     80976         0.99
## S_350 MSQ6_AG_381  52750    52500        1.00     51830     51938         0.99
## S_351 MSQ6_AG_382  73446    73177        1.00     72467     72411         0.99
## S_352 MSQ6_AG_383  69764    69491        1.00     68538     68691         0.99
## S_353 MSQ6_AG_384  61764    61575        1.00     60824     61048         0.99
## S_354 MSQ6_AG_385  52313    52053        1.00     51372     51569         0.99
## S_355 MSQ6_AG_386  51593    51313        0.99     50391     50675         0.98
## S_356 MSQ6_AG_387  46806    46542        0.99     45968     45972         0.99
## S_438 MSQ6_AG_490  66744    66455        1.00     65323     65479         0.98
## S_439 MSQ6_AG_491  65688    65409        1.00     64553     64591         0.99
## S_440 MSQ6_AG_492  42282    42132        1.00     41655     41507         0.99
## S_441 MSQ6_AG_493  70679    70422        1.00     69786     69632         0.99
## S_442 MSQ6_AG_494  64149    63878        1.00     62852     63031         0.98
## S_443 MSQ6_AG_495  66106    65772        0.99     65139     65014         0.99
## S_444 MSQ6_AG_496  67344    67025        1.00     65804     66358         0.98
## S_445 MSQ6_AG_497  47094    46707        0.99     45473     46244         0.97
## S_446 MSQ6_AG_498  67572    67267        1.00     65981     66428         0.98
## S_447 MSQ6_AG_499  67255    66994        1.00     65937     66416         0.98
## S_448 MSQ6_AG_500  57626    57445        1.00     56882     56862         0.99
## S_449 MSQ6_AG_501  69800    69629        1.00     69042     68944         0.99
## S_450 MSQ6_AG_502  66731    66475        1.00     65494     65868         0.99
## S_465 MSQ6_AG_517  69793    69596        1.00     68933     69150         0.99
## S_466 MSQ6_AG_518  73792    73415        0.99     72195     72677         0.98
## S_467 MSQ6_AG_519  64707    64401        1.00     63506     63635         0.99
## S_468 MSQ6_AG_520  69869    69411        0.99     67836     68545         0.98
## S_469 MSQ6_AG_521  22876    22575        0.99     21552     22319         0.95
## S_470 MSQ6_AG_522  70853    70436        0.99     69010     69691         0.98
## S_471 MSQ6_AG_523  66286    66030        1.00     64938     65385         0.98
## S_472 MSQ6_AG_524  63794    63539        1.00     62844     62818         0.99
## S_473 MSQ6_AG_525  76951    76627        1.00     76016     76028         0.99
## S_474 MSQ6_AG_526  61480    61249        1.00     60023     60475         0.98
## S_475 MSQ6_AG_527  65287    65039        1.00     64261     64598         0.99
## S_476 MSQ6_AG_528  62397    61982        0.99     60877     61304         0.98
## S_477 MSQ6_AG_529  40902    40453        0.99     39440     40103         0.97
## S_478 MSQ6_AG_530  67854    67566        1.00     66284     66723         0.98
##       denoisedR_pc merged merged_pc tabled chimera_out length_filtered
## S_219         0.99  37694      0.88  37694       26929           26929
## S_220         0.98  63881      0.90  63881       44125           44125
## S_221         0.99  67904      0.91  67904       49814           49814
## S_222         0.98  54905      0.89  54905       42532           42532
## S_223         0.99  68254      0.88  68254       49859           49859
## S_224         0.99  60502      0.85  60502       44187           44187
## S_225         0.99 173270      0.91 173270      124781          124781
## S_226         0.99  53780      0.90  53780       39110           39110
## S_227         0.99  65965      0.88  65965       47209           47209
## S_228         0.98  52534      0.91  52534       39238           39238
## S_229         0.99  59729      0.92  59729       49569           49569
## S_230         0.99  45226      0.90  45226       32999           32999
## S_231         0.99 113946      0.93 113946       89836           89836
## S_232         0.99  53809      0.90  53809       37806           37806
## S_233         0.99  58534      0.89  58534       44726           44726
## S_234         0.99  52635      0.91  52635       37761           37761
## S_251         0.99  42248      0.94  42248       36121           36121
## S_252         0.99  37345      0.89  37345       27664           27664
## S_253         0.99  55612      0.92  55612       39123           39123
## S_254         0.99  47517      0.93  47517       34414           34414
## S_255         0.99  55956      0.93  55956       41213           41213
## S_256         0.99  53508      0.93  53508       39884           39884
## S_257         0.99  48415      0.92  48415       35603           35603
## S_259         0.99  29649      0.88  29649       22393           22393
## S_260         0.99  53165      0.91  53165       37401           37401
## S_261         0.99  63752      0.91  63752       44075           44075
## S_262         0.99  42784      0.90  42784       33102           33102
## S_263         0.99  53184      0.90  53184       39177           39177
## S_264         0.99  51337      0.92  51337       36374           36371
## S_265         0.99  47000      0.91  47000       33074           33074
## S_267         0.98  33276      0.88  33276       25918           25918
## S_268         0.98  58828      0.90  58828       44286           44286
## S_269         0.99  55474      0.92  55474       41040           41040
## S_270         0.99  47109      0.93  47109       33787           33787
## S_271         0.99  62749      0.94  62749       49311           49311
## S_272         0.99  63342      0.93  63342       47625           47625
## S_273         0.99  48976      0.94  48976       38035           38035
## S_274         0.99  80614      0.93  80614       59937           59937
## S_275         0.99  42485      0.92  42485       30890           30890
## S_276         0.99  47698      0.92  47698       36986           36986
## S_277         0.98  52627      0.92  52627       40691           40691
## S_278         0.98  41374      0.91  41374       31897           31897
## S_279         0.99  64030      0.93  64030       47728           47728
## S_280         0.99  47420      0.93  47420       37306           37306
## S_281         0.99  43385      0.93  43385       32803           32803
## S_298         0.99  41344      0.93  41344       32872           32872
## S_299         0.99  32859      0.91  32859       25123           25123
## S_300         0.99  44429      0.93  44429       34215           34215
## S_301         0.99  50825      0.93  50825       39342           39342
## S_302         0.99  36206      0.94  36206       27584           27584
## S_303         0.99  51798      0.94  51798       36753           36753
## S_304         0.99  50705      0.93  50705       34090           34090
## S_305         0.99  43798      0.91  43798       32369           32369
## S_306         0.99  46319      0.92  46319       33486           33486
## S_308         0.98  51026      0.85  51026       37026           37026
## S_309         0.99  43102      0.89  43102       30284           30284
## S_310         0.99  33431      0.86  33431       24770           24770
## S_311         0.99  58799      0.88  58799       41706           41706
## S_312         0.99  47333      0.88  47333       34794           34794
## S_313         0.99  48093      0.94  48093       38160           38160
## S_315         0.99  50234      0.91  50234       37273           37273
## S_316         0.99  48036      0.93  48036       35747           35747
## S_317         0.99  50207      0.92  50207       36364           36364
## S_318         0.99  43142      0.93  43142       32689           32689
## S_319         0.99  46711      0.93  46711       35016           35016
## S_320         0.99  45769      0.92  45769       33891           33891
## S_321         0.99  48984      0.92  48984       35926           35926
## S_322         0.99  47055      0.94  47055       33868           33868
## S_323         0.99  50335      0.92  50335       38224           38224
## S_324         0.99  49621      0.95  49621       41368           41368
## S_325         0.99  39705      0.93  39705       29213           29213
## S_326         0.99  49886      0.92  49886       37943           37943
## S_327         0.99  41415      0.94  41415       33353           33353
## S_341         0.99  47667      0.94  47667       37388           37388
## S_342         0.99  36114      0.95  36114       30253           30253
## S_343         0.99  28431      0.92  28431       22578           22578
## S_344         0.99  52379      0.93  52379       41121           41121
## S_345         0.99  23818      0.90  23818       18120           18120
## S_346         0.99  38771      0.92  38771       30174           30174
## S_347         0.99  45983      0.94  45983       34059           34059
## S_348         0.99  46988      0.96  46988       36536           36536
## S_349         0.99  75663      0.94  75663       54422           54422
## S_350         0.99  48047      0.93  48047       33745           33745
## S_351         0.99  67645      0.93  67645       47662           47662
## S_352         0.99  63239      0.92  63239       45958           45958
## S_353         0.99  56477      0.93  56477       45184           45184
## S_354         0.99  48051      0.94  48051       34946           34946
## S_355         0.99  46169      0.92  46169       32861           32861
## S_356         0.99  41514      0.90  41514       29841           29841
## S_438         0.99  60561      0.93  60561       52584           52584
## S_439         0.99  60150      0.93  60150       50751           50751
## S_440         0.99  38560      0.93  38560       33283           33283
## S_441         0.99  66477      0.95  66477       58699           58699
## S_442         0.99  60377      0.96  60377       58881           58881
## S_443         0.99  61479      0.94  61479       50758           50758
## S_444         0.99  63124      0.96  63124       60863           60863
## S_445         0.99  42723      0.94  42723       36375           36375
## S_446         0.99  61694      0.94  61694       50918           50918
## S_447         0.99  63155      0.96  63155       58042           58042
## S_448         0.99  54043      0.95  54043       46975           46975
## S_449         0.99  65940      0.96  65940       57531           57531
## S_450         0.99  62165      0.95  62165       53894           53894
## S_465         0.99  67008      0.97  67008       63618           63618
## S_466         0.99  68679      0.95  68679       60000           60000
## S_467         0.99  59847      0.94  59847       49010           49010
## S_468         0.99  63359      0.93  63359       53096           53096
## S_469         0.99  19263      0.89  19263       15459           15459
## S_470         0.99  65289      0.95  65289       57063           57063
## S_471         0.99  58034      0.89  58034       43164           43164
## S_472         0.99  57373      0.91  57373       43837           43837
## S_473         0.99  72195      0.95  72195       56134           56134
## S_474         0.99  52673      0.88  52673       37386           37384
## S_475         0.99  60887      0.95  60887       48463           48463
## S_476         0.99  55620      0.91  55620       43056           43056
## S_477         0.99  35722      0.91  35722       27699           27699
## S_478         0.99  57328      0.86  57328       43663           43663
##       tabled_pc chimera_out_pc length_filtered_pc Sample reads remove
## S_219         1           0.71                  1  S_219 18616      s
## S_220         1           0.69                  1  S_220 27311      s
## S_221         1           0.73                  1  S_221 28127      s
## S_222         1           0.77                  1  S_222 23855      s
## S_223         1           0.73                  1  S_223 40793      s
## S_224         1           0.73                  1  S_224 37037      s
## S_225         1           0.72                  1  S_225 84206      s
## S_226         1           0.73                  1  S_226 24313      s
## S_227         1           0.72                  1  S_227 30945      s
## S_228         1           0.75                  1  S_228 17347      s
## S_229         1           0.83                  1  S_229 41233      s
## S_230         1           0.73                  1  S_230 17757      s
## S_231         1           0.79                  1  S_231 76765      s
## S_232         1           0.70                  1  S_232 24326      s
## S_233         1           0.76                  1  S_233 30767      s
## S_234         1           0.72                  1  S_234 23093      s
## S_251         1           0.85                  1  S_251 31371      s
## S_252         1           0.74                  1  S_252 19402      s
## S_253         1           0.70                  1  S_253 26676      s
## S_254         1           0.72                  1  S_254 23704      s
## S_255         1           0.74                  1  S_255 29094      s
## S_256         1           0.75                  1  S_256 27534      s
## S_257         1           0.74                  1  S_257 22057      s
## S_259         1           0.76                  1  S_259 15048      s
## S_260         1           0.70                  1  S_260 25086      s
## S_261         1           0.69                  1  S_261 30359      s
## S_262         1           0.77                  1  S_262 22352      s
## S_263         1           0.74                  1  S_263 34683      s
## S_264         1           0.71                  1  S_264 22771      s
## S_265         1           0.70                  1  S_265 26224      s
## S_267         1           0.78                  1  S_267 19363      s
## S_268         1           0.75                  1  S_268 33058      s
## S_269         1           0.74                  1  S_269 30595      s
## S_270         1           0.72                  1  S_270 25220      s
## S_271         1           0.79                  1  S_271 39359      s
## S_272         1           0.75                  1  S_272 34494      s
## S_273         1           0.78                  1  S_273 29387      s
## S_274         1           0.74                  1  S_274 48529      s
## S_275         1           0.73                  1  S_275 20899      s
## S_276         1           0.78                  1  S_276 26411      s
## S_277         1           0.77                  1  S_277 33269      s
## S_278         1           0.77                  1  S_278 24903      s
## S_279         1           0.75                  1  S_279 42828      s
## S_280         1           0.79                  1  S_280 24629      s
## S_281         1           0.76                  1  S_281 28530      s
## S_298         1           0.80                  1  S_298 23689      s
## S_299         1           0.76                  1  S_299 18327      s
## S_300         1           0.77                  1  S_300 23353      s
## S_301         1           0.77                  1  S_301 27068      s
## S_302         1           0.76                  1  S_302 16367      s
## S_303         1           0.71                  1  S_303 34241      s
## S_304         1           0.67                  1  S_304 32389      s
## S_305         1           0.74                  1  S_305 23822      s
## S_306         1           0.72                  1  S_306 20350      s
## S_308         1           0.73                  1  S_308 31462      s
## S_309         1           0.70                  1  S_309 14098      s
## S_310         1           0.74                  1  S_310 20222      s
## S_311         1           0.71                  1  S_311 32489      s
## S_312         1           0.74                  1  S_312 21782      s
## S_313         1           0.79                  1  S_313 24703      s
## S_315         1           0.74                  1  S_315 29733      s
## S_316         1           0.74                  1  S_316 25822      s
## S_317         1           0.72                  1  S_317 25042      s
## S_318         1           0.76                  1  S_318 22374      s
## S_319         1           0.75                  1  S_319 21462      s
## S_320         1           0.74                  1  S_320 24471      s
## S_321         1           0.73                  1  S_321 25050      s
## S_322         1           0.72                  1  S_322 16390      s
## S_323         1           0.76                  1  S_323 26938      s
## S_324         1           0.83                  1  S_324 24248      s
## S_325         1           0.74                  1  S_325 17733      s
## S_326         1           0.76                  1  S_326 34400      s
## S_327         1           0.81                  1  S_327 23792      s
## S_341         1           0.78                  1  S_341 26245      s
## S_342         1           0.84                  1  S_342 23351      s
## S_343         1           0.79                  1  S_343 16927      s
## S_344         1           0.79                  1  S_344 36013      s
## S_345         1           0.76                  1  S_345 16695      s
## S_346         1           0.78                  1  S_346 25765      s
## S_347         1           0.74                  1  S_347 24392      s
## S_348         1           0.78                  1  S_348 35828      s
## S_349         1           0.72                  1  S_349 42252      s
## S_350         1           0.70                  1  S_350 22422      s
## S_351         1           0.70                  1  S_351 35537      s
## S_352         1           0.73                  1  S_352 40142      s
## S_353         1           0.80                  1  S_353 39030      s
## S_354         1           0.73                  1  S_354 26885      s
## S_355         1           0.71                  1  S_355 23918      s
## S_356         1           0.72                  1  S_356 17693      s
## S_438         1           0.87                  1  S_438 42522      s
## S_439         1           0.84                  1  S_439 39448      s
## S_440         1           0.86                  1  S_440 25183      s
## S_441         1           0.88                  1  S_441 41289      s
## S_442         1           0.98                  1  S_442 45501      s
## S_443         1           0.83                  1  S_443 32594      s
## S_444         1           0.96                  1  S_444 47040      s
## S_445         1           0.85                  1  S_445 28263      s
## S_446         1           0.83                  1  S_446 39328      s
## S_447         1           0.92                  1  S_447 48814      s
## S_448         1           0.87                  1  S_448 34845      s
## S_449         1           0.87                  1  S_449 53167      s
## S_450         1           0.87                  1  S_450 46455      s
## S_465         1           0.95                  1  S_465 54378      s
## S_466         1           0.87                  1  S_466 46183      s
## S_467         1           0.82                  1  S_467 40486      s
## S_468         1           0.84                  1  S_468 44884      s
## S_469         1           0.80                  1  S_469 12799      s
## S_470         1           0.87                  1  S_470 47385      s
## S_471         1           0.74                  1  S_471 38843      s
## S_472         1           0.76                  1  S_472 30577      s
## S_473         1           0.78                  1  S_473 33737      s
## S_474         1           0.71                  1  S_474 36823      s
## S_475         1           0.80                  1  S_475 31171      s
## S_476         1           0.77                  1  S_476 32578      s
## S_477         1           0.78                  1  S_477 17380      s
## S_478         1           0.76                  1  S_478 37773      s
##       responsible_seq owner region_1  run experiment Fermentation facility
## S_219        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_220        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_221        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_222        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_223        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_224        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_225        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_226        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_227        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_228        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_229        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_230        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_231        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_232        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_233        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_234        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_251        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_252        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_253        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_254        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_255        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_256        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_257        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_259        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_260        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_261        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_262        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_263        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_264        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_265        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_267        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_268        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_269        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_270        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_271        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_272        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_273        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_274        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_275        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_276        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_277        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_278        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_279        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_280        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_281        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_298        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_299        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_300        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_301        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_302        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_303        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_304        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_305        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_306        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_308        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_309        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_310        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_311        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_312        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_313        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_315        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_316        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_317        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_318        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_319        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_320        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_321        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_322        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_323        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_324        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_325        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_326        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_327        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_341        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_342        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_343        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_344        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_345        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_346        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_347        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_348        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_349        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_350        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_351        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_352        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_353        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_354        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_355        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_356        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_438        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_439        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_440        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_441        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_442        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_443        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_444        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_445        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_446        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_447        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_448        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_449        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_450        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_465        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_466        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_467        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_468        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_469        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_470        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_471        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_472        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_473        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_474        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_475        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_476        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_477        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
## S_478        Annelies    AG       V4 MSQ6     mIMT_2         <NA>   Scharl
##          day sample_type reactor day_IR day_fermentation treatment_period
## S_219 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_220 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_221 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_222 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_223 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_224 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_225 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_226 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_227 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_228 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_229 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_230 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_231 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_232 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_233 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_234 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_251 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_252 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_253 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_254 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_255 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_256 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_257 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_259 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_260 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_261 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_262 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_263 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_264 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_265 Day_-1       Feces    <NA>     NA             <NA>             <NA>
## S_267 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_268 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_269 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_270 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_271 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_272 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_273 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_274 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_275 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_276 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_277 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_278 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_279 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_280 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_281 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_298 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_299 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_300 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_301 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_302 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_303 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_304 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_305 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_306 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_308 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_309 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_310 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_311 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_312 Day_12       Feces    <NA>     NA             <NA>             <NA>
## S_313 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_315 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_316 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_317 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_318 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_319 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_320 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_321 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_322 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_323 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_324 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_325 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_326 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_327 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_341 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_342 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_343 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_344 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_345 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_346 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_347 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_348 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_349 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_350 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_351 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_352 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_353 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_354 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_355 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_356 Day_14       Feces    <NA>     NA             <NA>             <NA>
## S_438  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_439  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_440  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_441  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_442  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_443  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_444  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_445  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_446  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_447  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_448  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_449  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_450  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_465  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_466  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_467  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_468  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_469  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_470  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_471  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_472  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_473  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_474  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_475  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_476  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_477  Day_6       Feces    <NA>     NA             <NA>             <NA>
## S_478  Day_6       Feces    <NA>     NA             <NA>             <NA>
##       sample_label treatment_exvivo treatment mouse_label    BW   BW_percent
## S_219         <NA>             <NA>  Eubiotic         m_8 20.55  0.000000000
## S_220         <NA>             <NA>  Eubiotic        m_19 21.32  0.000000000
## S_221         <NA>             <NA>  Eubiotic        m_34 20.11  0.000000000
## S_222         <NA>             <NA>  Eubiotic        m_26 21.90  0.000000000
## S_223         <NA>             <NA>  Eubiotic        m_39 20.94  0.000000000
## S_224         <NA>             <NA>  Eubiotic        m_46 21.12  0.000000000
## S_225         <NA>             <NA>  Eubiotic        m_48 20.35  0.000000000
## S_226         <NA>             <NA>  Eubiotic        m_56 20.51  0.000000000
## S_227         <NA>             <NA> Dysbiotic        m_29 20.95  0.000000000
## S_228         <NA>             <NA> Dysbiotic        m_13 20.35  0.000000000
## S_229         <NA>             <NA> Dysbiotic        m_18 21.55  0.000000000
## S_230         <NA>             <NA> Dysbiotic        m_27 23.18  0.000000000
## S_231         <NA>             <NA> Dysbiotic        m_32 20.78  0.000000000
## S_232         <NA>             <NA> Dysbiotic        m_49 20.90  0.000000000
## S_233         <NA>             <NA> Dysbiotic        m_54 19.70  0.000000000
## S_234         <NA>             <NA> Dysbiotic        m_35 20.53  0.000000000
## S_251         <NA>             <NA>       DSS         m_7 20.98  0.000000000
## S_252         <NA>             <NA>       DSS        m_10 20.48  0.000000000
## S_253         <NA>             <NA>       DSS        m_11 20.47  0.000000000
## S_254         <NA>             <NA>       DSS        m_24 22.91  0.000000000
## S_255         <NA>             <NA>       DSS        m_28 21.30  0.000000000
## S_256         <NA>             <NA>       DSS        m_36 22.47  0.000000000
## S_257         <NA>             <NA>       DSS        m_51 21.92  0.000000000
## S_259         <NA>             <NA>       H2O         m_2 21.09  0.000000000
## S_260         <NA>             <NA>       H2O        m_21 21.31  0.000000000
## S_261         <NA>             <NA>       H2O        m_22 20.45  0.000000000
## S_262         <NA>             <NA>       H2O        m_25 22.00  0.000000000
## S_263         <NA>             <NA>       H2O        m_37 20.38  0.000000000
## S_264         <NA>             <NA>       H2O        m_44 21.60  0.000000000
## S_265         <NA>             <NA>       H2O        m_52 21.75  0.000000000
## S_267         <NA>             <NA>  Eubiotic         m_8 18.24  0.112408759
## S_268         <NA>             <NA>  Eubiotic        m_19 18.64  0.125703565
## S_269         <NA>             <NA>  Eubiotic        m_34 19.75  0.017901542
## S_270         <NA>             <NA>  Eubiotic        m_26 22.45 -0.025114155
## S_271         <NA>             <NA>  Eubiotic        m_39 20.86  0.003820439
## S_272         <NA>             <NA>  Eubiotic        m_46 20.12  0.047348485
## S_273         <NA>             <NA>  Eubiotic        m_48 19.75  0.029484029
## S_274         <NA>             <NA> Dysbiotic        m_29 17.62  0.158949881
## S_275         <NA>             <NA> Dysbiotic        m_13 17.60  0.135135135
## S_276         <NA>             <NA> Dysbiotic        m_18 17.57  0.184686775
## S_277         <NA>             <NA> Dysbiotic        m_27 20.47  0.116911130
## S_278         <NA>             <NA> Dysbiotic        m_32 17.32  0.166506256
## S_279         <NA>             <NA> Dysbiotic        m_49 17.95  0.141148325
## S_280         <NA>             <NA> Dysbiotic        m_54 17.25  0.124365482
## S_281         <NA>             <NA> Dysbiotic        m_35 18.45  0.101315149
## S_298         <NA>             <NA>       DSS         m_7 18.25  0.130123928
## S_299         <NA>             <NA>       DSS        m_10 17.35  0.152832031
## S_300         <NA>             <NA>       DSS        m_11 16.85  0.176844162
## S_301         <NA>             <NA>       DSS        m_24 20.01  0.126582278
## S_302         <NA>             <NA>       DSS        m_28 18.12  0.149295775
## S_303         <NA>             <NA>       DSS        m_36 17.07  0.240320427
## S_304         <NA>             <NA>       DSS        m_51 16.07  0.266879562
## S_305         <NA>             <NA>       H2O         m_2 20.70 -0.018492176
## S_306         <NA>             <NA>       H2O        m_21 21.54  0.010793055
## S_308         <NA>             <NA>       H2O        m_25 21.86 -0.006363636
## S_309         <NA>             <NA>       H2O        m_37 21.16  0.038272816
## S_310         <NA>             <NA>       H2O        m_44 22.97  0.063425926
## S_311         <NA>             <NA>       H2O        m_52 22.07  0.014712644
## S_312         <NA>             <NA>       H2O        m_53 21.47 -0.029823769
## S_313         <NA>             <NA>  Eubiotic        m_19 21.13  0.008911820
## S_315         <NA>             <NA>  Eubiotic        m_34 21.03 -0.045748384
## S_316         <NA>             <NA>  Eubiotic        m_26 21.34  0.025570776
## S_317         <NA>             <NA>  Eubiotic        m_39 22.88 -0.092645654
## S_318         <NA>             <NA>  Eubiotic        m_46 21.81 -0.032670455
## S_319         <NA>             <NA>  Eubiotic        m_48 21.31 -0.047174447
## S_320         <NA>             <NA> Dysbiotic        m_29 18.02  0.139856802
## S_321         <NA>             <NA> Dysbiotic        m_13 18.23  0.104176904
## S_322         <NA>             <NA> Dysbiotic        m_18 18.01  0.164269142
## S_323         <NA>             <NA> Dysbiotic        m_27 20.95  0.096203624
## S_324         <NA>             <NA> Dysbiotic        m_32 17.85  0.141000962
## S_325         <NA>             <NA> Dysbiotic        m_49 18.35  0.122009569
## S_326         <NA>             <NA> Dysbiotic        m_54 18.52  0.059898477
## S_327         <NA>             <NA> Dysbiotic        m_35 18.95  0.076960546
## S_341         <NA>             <NA>       DSS         m_7 18.21  0.132030505
## S_342         <NA>             <NA>       DSS        m_10 17.58  0.141601563
## S_343         <NA>             <NA>       DSS        m_11 17.65  0.137762579
## S_344         <NA>             <NA>       DSS        m_24 20.23  0.116979485
## S_345         <NA>             <NA>       DSS        m_28 18.96  0.109859155
## S_346         <NA>             <NA>       DSS        m_36 18.02  0.198041834
## S_347         <NA>             <NA>       DSS        m_51 16.96  0.226277372
## S_348         <NA>             <NA>       DSS        m_55 18.02  0.121404193
## S_349         <NA>             <NA>       H2O         m_2 21.28  0.009009009
## S_350         <NA>             <NA>       H2O        m_21 21.72  0.019239794
## S_351         <NA>             <NA>       H2O        m_22 21.70  0.061124694
## S_352         <NA>             <NA>       H2O        m_25 21.79 -0.009545455
## S_353         <NA>             <NA>       H2O        m_37 21.72  0.065750736
## S_354         <NA>             <NA>       H2O        m_44 23.77  0.100462963
## S_355         <NA>             <NA>       H2O        m_52 21.75  0.000000000
## S_356         <NA>             <NA>       H2O        m_53 22.26  0.005874379
## S_438         <NA>             <NA>  Eubiotic        m_19 20.48  0.039399625
## S_439         <NA>             <NA>  Eubiotic        m_26 20.35  0.070776256
## S_440         <NA>             <NA>  Eubiotic        m_46 20.77  0.016571970
## S_441         <NA>             <NA>  Eubiotic        m_48 19.90  0.022113022
## S_442         <NA>             <NA>  Eubiotic        m_56 19.50  0.049244271
## S_443         <NA>             <NA> Dysbiotic        m_29 21.32 -0.017661098
## S_444         <NA>             <NA> Dysbiotic        m_13 20.59 -0.011793612
## S_445         <NA>             <NA> Dysbiotic        m_18 21.06  0.022737819
## S_446         <NA>             <NA> Dysbiotic        m_27 22.95  0.009922347
## S_447         <NA>             <NA> Dysbiotic        m_32 20.91 -0.006256015
## S_448         <NA>             <NA> Dysbiotic        m_49 20.41  0.023444976
## S_449         <NA>             <NA> Dysbiotic        m_54 19.93 -0.011675127
## S_450         <NA>             <NA> Dysbiotic        m_35 20.81 -0.013638578
## S_465         <NA>             <NA>       DSS         m_7 20.98  0.000000000
## S_466         <NA>             <NA>       DSS        m_10 20.42  0.002929687
## S_467         <NA>             <NA>       DSS        m_11 19.93  0.026380068
## S_468         <NA>             <NA>       DSS        m_24 23.81 -0.039284155
## S_469         <NA>             <NA>       DSS        m_28 20.73  0.026760563
## S_470         <NA>             <NA>       DSS        m_36 21.79  0.030262572
## S_471         <NA>             <NA>       H2O         m_2 20.65 -0.020862968
## S_472         <NA>             <NA>       H2O        m_21 21.20 -0.005161896
## S_473         <NA>             <NA>       H2O        m_22 20.90  0.022004890
## S_474         <NA>             <NA>       H2O        m_25 21.68 -0.014545455
## S_475         <NA>             <NA>       H2O        m_37 21.68  0.063788027
## S_476         <NA>             <NA>       H2O        m_44 22.45  0.039351852
## S_477         <NA>             <NA>       H2O        m_52 20.76 -0.045517241
## S_478         <NA>             <NA>       H2O        m_53 21.74 -0.017623136
##       BW_delta DAI colon_lenght spleen_weight histo_colo_infilt
## S_219     0.00   0           NA            NA                NA
## S_220     0.00   0           NA            NA                NA
## S_221     0.00   0           NA            NA                NA
## S_222     0.00   0           NA            NA                NA
## S_223     0.00   0           NA            NA                NA
## S_224     0.00   0           NA            NA                NA
## S_225     0.00   0           NA            NA                NA
## S_226     0.00   0           NA            NA                NA
## S_227     0.00   0           NA            NA                NA
## S_228     0.00   0           NA            NA                NA
## S_229     0.00   0           NA            NA                NA
## S_230     0.00   0           NA            NA                NA
## S_231     0.00   0           NA            NA                NA
## S_232     0.00   0           NA            NA                NA
## S_233     0.00   0           NA            NA                NA
## S_234     0.00   0           NA            NA                NA
## S_251     0.00   0           NA            NA                NA
## S_252     0.00   0           NA            NA                NA
## S_253     0.00   0           NA            NA                NA
## S_254     0.00   0           NA            NA                NA
## S_255     0.00   0           NA            NA                NA
## S_256     0.00   0           NA            NA                NA
## S_257     0.00   0           NA            NA                NA
## S_259     0.00   0           NA            NA                NA
## S_260     0.00   0           NA            NA                NA
## S_261     0.00   0           NA            NA                NA
## S_262     0.00   0           NA            NA                NA
## S_263     0.00   0           NA            NA                NA
## S_264     0.00   0           NA            NA                NA
## S_265     0.00   0           NA            NA                NA
## S_267     0.20   1           NA            NA                NA
## S_268     0.64   1           NA            NA                NA
## S_269     0.90   0           NA            NA                NA
## S_270    -0.15   1           NA            NA                NA
## S_271     1.81   0           NA            NA                NA
## S_272     0.70   0           NA            NA                NA
## S_273     0.71   0           NA            NA                NA
## S_274     0.26   3           NA            NA                NA
## S_275     0.05   3           NA            NA                NA
## S_276     0.39   3           NA            NA                NA
## S_277    -0.46   1           NA            NA                NA
## S_278     0.24   3           NA            NA                NA
## S_279     0.18   3           NA            NA                NA
## S_280     0.20   3           NA            NA                NA
## S_281    -0.10   3           NA            NA                NA
## S_298     0.30   3           NA            NA                NA
## S_299     0.25   3           NA            NA                NA
## S_300    -0.01   4           NA            NA                NA
## S_301    -0.35   2           NA            NA                NA
## S_302    -0.05   3           NA            NA                NA
## S_303    -0.19   4           NA            NA                NA
## S_304     0.20   4           NA            NA                NA
## S_305    -0.02   0           NA            NA                NA
## S_306     0.74   0           NA            NA                NA
## S_308    -1.88   0           NA            NA                NA
## S_309    -0.22   0           NA            NA                NA
## S_310     0.12   0           NA            NA                NA
## S_311     0.67   0           NA            NA                NA
## S_312    -0.87   0           NA            NA                NA
## S_313     1.32   0          5.2         0.105               1.0
## S_315     0.84   0          5.4         0.086               2.5
## S_316    -0.85   0          6.0         0.095               2.5
## S_317     2.42   0          5.2         0.101               1.0
## S_318     0.59   0          6.1         0.128               1.5
## S_319     0.49   0          5.8         0.101               2.0
## S_320     0.17   2          4.0         0.154               3.0
## S_321     0.38   2          4.2         0.162               4.0
## S_322     0.45   3          4.5         0.144               3.0
## S_323     0.94   0          4.2         0.128               3.5
## S_324     0.50   2          4.6         0.152               2.5
## S_325     0.23   2          4.9         0.129               3.0
## S_326     0.52   0          5.1         0.136               3.0
## S_327     0.41   0          4.3         0.148               1.5
## S_341     0.07   2          4.2         0.147               4.0
## S_342    -0.07   3          4.7         0.123               3.5
## S_343     0.09   2          4.3         0.125               3.0
## S_344    -0.40   3          4.6         0.142               2.5
## S_345     0.42   3          4.4         0.110               4.0
## S_346     0.80   3          4.9         0.183               3.5
## S_347     0.72   4          4.3         0.165               3.5
## S_348     0.06   2          4.2         0.145               4.0
## S_349     0.58   0          7.8         0.065               0.5
## S_350     0.18   0          6.2         0.093               0.0
## S_351     0.54   0          6.9         0.074               1.0
## S_352    -0.07   0          7.1         0.071               0.0
## S_353     0.56   0          6.8         0.067               0.0
## S_354     0.80   0          6.7         0.082               0.0
## S_355    -0.26   0          6.5         0.065               0.0
## S_356     0.79   0          6.9         0.065               0.0
## S_438    -0.93   1           NA            NA                NA
## S_439    -1.61   2           NA            NA                NA
## S_440    -0.49   1           NA            NA                NA
## S_441    -0.43   1           NA            NA                NA
## S_442    -1.12   1           NA            NA                NA
## S_443     0.39   2           NA            NA                NA
## S_444    -0.79   1           NA            NA                NA
## S_445    -1.04   1           NA            NA                NA
## S_446    -0.21   0           NA            NA                NA
## S_447    -0.15   1           NA            NA                NA
## S_448    -1.00   2           NA            NA                NA
## S_449    -0.86   1           NA            NA                NA
## S_450    -0.60   1           NA            NA                NA
## S_465    -0.93   1           NA            NA                NA
## S_466    -0.69   1           NA            NA                NA
## S_467    -0.87   1           NA            NA                NA
## S_468     0.46   0           NA            NA                NA
## S_469    -0.70   1           NA            NA                NA
## S_470    -0.33   2           NA            NA                NA
## S_471    -0.34   0           NA            NA                NA
## S_472     0.25   0           NA            NA                NA
## S_473    -1.05   0           NA            NA                NA
## S_474     0.38   0           NA            NA                NA
## S_475     0.43   0           NA            NA                NA
## S_476     0.52   0           NA            NA                NA
## S_477    -0.19   0           NA            NA                NA
## S_478     0.18   0           NA            NA                NA
##       histo_colo_damage histo_colo_total infiltr_CD3 infiltr_Gr1 CM_SCFA
## S_219                NA               NA          NA          NA      NA
## S_220                NA               NA          NA          NA      NA
## S_221                NA               NA          NA          NA      NA
## S_222                NA               NA          NA          NA      NA
## S_223                NA               NA          NA          NA      NA
## S_224                NA               NA          NA          NA      NA
## S_225                NA               NA          NA          NA      NA
## S_226                NA               NA          NA          NA      NA
## S_227                NA               NA          NA          NA      NA
## S_228                NA               NA          NA          NA      NA
## S_229                NA               NA          NA          NA      NA
## S_230                NA               NA          NA          NA      NA
## S_231                NA               NA          NA          NA      NA
## S_232                NA               NA          NA          NA      NA
## S_233                NA               NA          NA          NA      NA
## S_234                NA               NA          NA          NA      NA
## S_251                NA               NA          NA          NA      NA
## S_252                NA               NA          NA          NA      NA
## S_253                NA               NA          NA          NA      NA
## S_254                NA               NA          NA          NA      NA
## S_255                NA               NA          NA          NA      NA
## S_256                NA               NA          NA          NA      NA
## S_257                NA               NA          NA          NA      NA
## S_259                NA               NA          NA          NA      NA
## S_260                NA               NA          NA          NA      NA
## S_261                NA               NA          NA          NA      NA
## S_262                NA               NA          NA          NA      NA
## S_263                NA               NA          NA          NA      NA
## S_264                NA               NA          NA          NA      NA
## S_265                NA               NA          NA          NA      NA
## S_267                NA               NA          NA          NA      NA
## S_268                NA               NA          NA          NA      NA
## S_269                NA               NA          NA          NA      NA
## S_270                NA               NA          NA          NA      NA
## S_271                NA               NA          NA          NA      NA
## S_272                NA               NA          NA          NA      NA
## S_273                NA               NA          NA          NA      NA
## S_274                NA               NA          NA          NA      NA
## S_275                NA               NA          NA          NA      NA
## S_276                NA               NA          NA          NA      NA
## S_277                NA               NA          NA          NA      NA
## S_278                NA               NA          NA          NA      NA
## S_279                NA               NA          NA          NA      NA
## S_280                NA               NA          NA          NA      NA
## S_281                NA               NA          NA          NA      NA
## S_298                NA               NA          NA          NA      NA
## S_299                NA               NA          NA          NA      NA
## S_300                NA               NA          NA          NA      NA
## S_301                NA               NA          NA          NA      NA
## S_302                NA               NA          NA          NA      NA
## S_303                NA               NA          NA          NA      NA
## S_304                NA               NA          NA          NA      NA
## S_305                NA               NA          NA          NA      NA
## S_306                NA               NA          NA          NA      NA
## S_308                NA               NA          NA          NA      NA
## S_309                NA               NA          NA          NA      NA
## S_310                NA               NA          NA          NA      NA
## S_311                NA               NA          NA          NA      NA
## S_312                NA               NA          NA          NA      NA
## S_313               1.5              2.5          50          65      NA
## S_315               2.0              4.5          49          27      NA
## S_316               2.0              4.5          26          30      NA
## S_317               1.0              2.0          72          12      NA
## S_318               1.0              2.5          26          76      NA
## S_319               2.0              4.0          46          10      NA
## S_320               3.5              6.5         125         130      NA
## S_321               4.0              8.0          89          76      NA
## S_322               2.5              5.5          74          80      NA
## S_323               3.0              6.5         134         145      NA
## S_324               2.5              5.0         186         155      NA
## S_325               4.0              7.0          92         101      NA
## S_326               3.0              6.0         119          89      NA
## S_327               2.0              3.5          78         123      NA
## S_341               4.0              8.0         135          98      NA
## S_342               3.5              7.0         100         142      NA
## S_343               4.0              7.0         144          79      NA
## S_344               3.0              5.5         132          85      NA
## S_345               4.0              8.0         112          97      NA
## S_346               4.0              7.5          97          79      NA
## S_347               4.0              7.5         122          87      NA
## S_348               3.0              7.0         154         115      NA
## S_349               0.0              0.5          15           5      NA
## S_350               0.0              0.0          17           8      NA
## S_351               0.5              1.5          32          14      NA
## S_352               0.0              0.0          18           2      NA
## S_353               0.0              0.0          25          17      NA
## S_354               0.0              0.0          10          33      NA
## S_355               0.0              0.0           8          10      NA
## S_356               0.5              0.5          14           6      NA
## S_438                NA               NA          NA          NA      NA
## S_439                NA               NA          NA          NA      NA
## S_440                NA               NA          NA          NA      NA
## S_441                NA               NA          NA          NA      NA
## S_442                NA               NA          NA          NA      NA
## S_443                NA               NA          NA          NA      NA
## S_444                NA               NA          NA          NA      NA
## S_445                NA               NA          NA          NA      NA
## S_446                NA               NA          NA          NA      NA
## S_447                NA               NA          NA          NA      NA
## S_448                NA               NA          NA          NA      NA
## S_449                NA               NA          NA          NA      NA
## S_450                NA               NA          NA          NA      NA
## S_465                NA               NA          NA          NA      NA
## S_466                NA               NA          NA          NA      NA
## S_467                NA               NA          NA          NA      NA
## S_468                NA               NA          NA          NA      NA
## S_469                NA               NA          NA          NA      NA
## S_470                NA               NA          NA          NA      NA
## S_471                NA               NA          NA          NA      NA
## S_472                NA               NA          NA          NA      NA
## S_473                NA               NA          NA          NA      NA
## S_474                NA               NA          NA          NA      NA
## S_475                NA               NA          NA          NA      NA
## S_476                NA               NA          NA          NA      NA
## S_477                NA               NA          NA          NA      NA
## S_478                NA               NA          NA          NA      NA
##       CM_BSCFA CM_total CM_acetate CM_propioNAte CM_butyrate CM_succiNAte
## S_219       NA       NA         NA            NA          NA           NA
## S_220       NA       NA         NA            NA          NA           NA
## S_221       NA       NA         NA            NA          NA           NA
## S_222       NA       NA         NA            NA          NA           NA
## S_223       NA       NA         NA            NA          NA           NA
## S_224       NA       NA         NA            NA          NA           NA
## S_225       NA       NA         NA            NA          NA           NA
## S_226       NA       NA         NA            NA          NA           NA
## S_227       NA       NA         NA            NA          NA           NA
## S_228       NA       NA         NA            NA          NA           NA
## S_229       NA       NA         NA            NA          NA           NA
## S_230       NA       NA         NA            NA          NA           NA
## S_231       NA       NA         NA            NA          NA           NA
## S_232       NA       NA         NA            NA          NA           NA
## S_233       NA       NA         NA            NA          NA           NA
## S_234       NA       NA         NA            NA          NA           NA
## S_251       NA       NA         NA            NA          NA           NA
## S_252       NA       NA         NA            NA          NA           NA
## S_253       NA       NA         NA            NA          NA           NA
## S_254       NA       NA         NA            NA          NA           NA
## S_255       NA       NA         NA            NA          NA           NA
## S_256       NA       NA         NA            NA          NA           NA
## S_257       NA       NA         NA            NA          NA           NA
## S_259       NA       NA         NA            NA          NA           NA
## S_260       NA       NA         NA            NA          NA           NA
## S_261       NA       NA         NA            NA          NA           NA
## S_262       NA       NA         NA            NA          NA           NA
## S_263       NA       NA         NA            NA          NA           NA
## S_264       NA       NA         NA            NA          NA           NA
## S_265       NA       NA         NA            NA          NA           NA
## S_267       NA       NA         NA            NA          NA           NA
## S_268       NA       NA         NA            NA          NA           NA
## S_269       NA       NA         NA            NA          NA           NA
## S_270       NA       NA         NA            NA          NA           NA
## S_271       NA       NA         NA            NA          NA           NA
## S_272       NA       NA         NA            NA          NA           NA
## S_273       NA       NA         NA            NA          NA           NA
## S_274       NA       NA         NA            NA          NA           NA
## S_275       NA       NA         NA            NA          NA           NA
## S_276       NA       NA         NA            NA          NA           NA
## S_277       NA       NA         NA            NA          NA           NA
## S_278       NA       NA         NA            NA          NA           NA
## S_279       NA       NA         NA            NA          NA           NA
## S_280       NA       NA         NA            NA          NA           NA
## S_281       NA       NA         NA            NA          NA           NA
## S_298       NA       NA         NA            NA          NA           NA
## S_299       NA       NA         NA            NA          NA           NA
## S_300       NA       NA         NA            NA          NA           NA
## S_301       NA       NA         NA            NA          NA           NA
## S_302       NA       NA         NA            NA          NA           NA
## S_303       NA       NA         NA            NA          NA           NA
## S_304       NA       NA         NA            NA          NA           NA
## S_305       NA       NA         NA            NA          NA           NA
## S_306       NA       NA         NA            NA          NA           NA
## S_308       NA       NA         NA            NA          NA           NA
## S_309       NA       NA         NA            NA          NA           NA
## S_310       NA       NA         NA            NA          NA           NA
## S_311       NA       NA         NA            NA          NA           NA
## S_312       NA       NA         NA            NA          NA           NA
## S_313       NA       NA         NA            NA          NA           NA
## S_315       NA       NA         NA            NA          NA           NA
## S_316       NA       NA         NA            NA          NA           NA
## S_317       NA       NA         NA            NA          NA           NA
## S_318       NA       NA         NA            NA          NA           NA
## S_319       NA       NA         NA            NA          NA           NA
## S_320       NA       NA         NA            NA          NA           NA
## S_321       NA       NA         NA            NA          NA           NA
## S_322       NA       NA         NA            NA          NA           NA
## S_323       NA       NA         NA            NA          NA           NA
## S_324       NA       NA         NA            NA          NA           NA
## S_325       NA       NA         NA            NA          NA           NA
## S_326       NA       NA         NA            NA          NA           NA
## S_327       NA       NA         NA            NA          NA           NA
## S_341       NA       NA         NA            NA          NA           NA
## S_342       NA       NA         NA            NA          NA           NA
## S_343       NA       NA         NA            NA          NA           NA
## S_344       NA       NA         NA            NA          NA           NA
## S_345       NA       NA         NA            NA          NA           NA
## S_346       NA       NA         NA            NA          NA           NA
## S_347       NA       NA         NA            NA          NA           NA
## S_348       NA       NA         NA            NA          NA           NA
## S_349       NA       NA         NA            NA          NA           NA
## S_350       NA       NA         NA            NA          NA           NA
## S_351       NA       NA         NA            NA          NA           NA
## S_352       NA       NA         NA            NA          NA           NA
## S_353       NA       NA         NA            NA          NA           NA
## S_354       NA       NA         NA            NA          NA           NA
## S_355       NA       NA         NA            NA          NA           NA
## S_356       NA       NA         NA            NA          NA           NA
## S_438       NA       NA         NA            NA          NA           NA
## S_439       NA       NA         NA            NA          NA           NA
## S_440       NA       NA         NA            NA          NA           NA
## S_441       NA       NA         NA            NA          NA           NA
## S_442       NA       NA         NA            NA          NA           NA
## S_443       NA       NA         NA            NA          NA           NA
## S_444       NA       NA         NA            NA          NA           NA
## S_445       NA       NA         NA            NA          NA           NA
## S_446       NA       NA         NA            NA          NA           NA
## S_447       NA       NA         NA            NA          NA           NA
## S_448       NA       NA         NA            NA          NA           NA
## S_449       NA       NA         NA            NA          NA           NA
## S_450       NA       NA         NA            NA          NA           NA
## S_465       NA       NA         NA            NA          NA           NA
## S_466       NA       NA         NA            NA          NA           NA
## S_467       NA       NA         NA            NA          NA           NA
## S_468       NA       NA         NA            NA          NA           NA
## S_469       NA       NA         NA            NA          NA           NA
## S_470       NA       NA         NA            NA          NA           NA
## S_471       NA       NA         NA            NA          NA           NA
## S_472       NA       NA         NA            NA          NA           NA
## S_473       NA       NA         NA            NA          NA           NA
## S_474       NA       NA         NA            NA          NA           NA
## S_475       NA       NA         NA            NA          NA           NA
## S_476       NA       NA         NA            NA          NA           NA
## S_477       NA       NA         NA            NA          NA           NA
## S_478       NA       NA         NA            NA          NA           NA
##       CM_lactate CB_total CB_intact CB_permeable QPCR_total QPCR_Akkermansia
## S_219         NA       NA        NA           NA         NA               NA
## S_220         NA       NA        NA           NA         NA               NA
## S_221         NA       NA        NA           NA         NA               NA
## S_222         NA       NA        NA           NA         NA               NA
## S_223         NA       NA        NA           NA         NA               NA
## S_224         NA       NA        NA           NA         NA               NA
## S_225         NA       NA        NA           NA         NA               NA
## S_226         NA       NA        NA           NA         NA               NA
## S_227         NA       NA        NA           NA         NA               NA
## S_228         NA       NA        NA           NA         NA               NA
## S_229         NA       NA        NA           NA         NA               NA
## S_230         NA       NA        NA           NA         NA               NA
## S_231         NA       NA        NA           NA         NA               NA
## S_232         NA       NA        NA           NA         NA               NA
## S_233         NA       NA        NA           NA         NA               NA
## S_234         NA       NA        NA           NA         NA               NA
## S_251         NA       NA        NA           NA         NA               NA
## S_252         NA       NA        NA           NA         NA               NA
## S_253         NA       NA        NA           NA         NA               NA
## S_254         NA       NA        NA           NA         NA               NA
## S_255         NA       NA        NA           NA         NA               NA
## S_256         NA       NA        NA           NA         NA               NA
## S_257         NA       NA        NA           NA         NA               NA
## S_259         NA       NA        NA           NA         NA               NA
## S_260         NA       NA        NA           NA         NA               NA
## S_261         NA       NA        NA           NA         NA               NA
## S_262         NA       NA        NA           NA         NA               NA
## S_263         NA       NA        NA           NA         NA               NA
## S_264         NA       NA        NA           NA         NA               NA
## S_265         NA       NA        NA           NA         NA               NA
## S_267         NA       NA        NA           NA         NA               NA
## S_268         NA       NA        NA           NA         NA               NA
## S_269         NA       NA        NA           NA         NA               NA
## S_270         NA       NA        NA           NA         NA               NA
## S_271         NA       NA        NA           NA         NA               NA
## S_272         NA       NA        NA           NA         NA               NA
## S_273         NA       NA        NA           NA         NA               NA
## S_274         NA       NA        NA           NA         NA               NA
## S_275         NA       NA        NA           NA         NA               NA
## S_276         NA       NA        NA           NA         NA               NA
## S_277         NA       NA        NA           NA         NA               NA
## S_278         NA       NA        NA           NA         NA               NA
## S_279         NA       NA        NA           NA         NA               NA
## S_280         NA       NA        NA           NA         NA               NA
## S_281         NA       NA        NA           NA         NA               NA
## S_298         NA       NA        NA           NA         NA               NA
## S_299         NA       NA        NA           NA         NA               NA
## S_300         NA       NA        NA           NA         NA               NA
## S_301         NA       NA        NA           NA         NA               NA
## S_302         NA       NA        NA           NA         NA               NA
## S_303         NA       NA        NA           NA         NA               NA
## S_304         NA       NA        NA           NA         NA               NA
## S_305         NA       NA        NA           NA         NA               NA
## S_306         NA       NA        NA           NA         NA               NA
## S_308         NA       NA        NA           NA         NA               NA
## S_309         NA       NA        NA           NA         NA               NA
## S_310         NA       NA        NA           NA         NA               NA
## S_311         NA       NA        NA           NA         NA               NA
## S_312         NA       NA        NA           NA         NA               NA
## S_313         NA       NA        NA           NA         NA               NA
## S_315         NA       NA        NA           NA         NA               NA
## S_316         NA       NA        NA           NA         NA               NA
## S_317         NA       NA        NA           NA         NA               NA
## S_318         NA       NA        NA           NA         NA               NA
## S_319         NA       NA        NA           NA         NA               NA
## S_320         NA       NA        NA           NA         NA               NA
## S_321         NA       NA        NA           NA         NA               NA
## S_322         NA       NA        NA           NA         NA               NA
## S_323         NA       NA        NA           NA         NA               NA
## S_324         NA       NA        NA           NA         NA               NA
## S_325         NA       NA        NA           NA         NA               NA
## S_326         NA       NA        NA           NA         NA               NA
## S_327         NA       NA        NA           NA         NA               NA
## S_341         NA       NA        NA           NA         NA               NA
## S_342         NA       NA        NA           NA         NA               NA
## S_343         NA       NA        NA           NA         NA               NA
## S_344         NA       NA        NA           NA         NA               NA
## S_345         NA       NA        NA           NA         NA               NA
## S_346         NA       NA        NA           NA         NA               NA
## S_347         NA       NA        NA           NA         NA               NA
## S_348         NA       NA        NA           NA         NA               NA
## S_349         NA       NA        NA           NA         NA               NA
## S_350         NA       NA        NA           NA         NA               NA
## S_351         NA       NA        NA           NA         NA               NA
## S_352         NA       NA        NA           NA         NA               NA
## S_353         NA       NA        NA           NA         NA               NA
## S_354         NA       NA        NA           NA         NA               NA
## S_355         NA       NA        NA           NA         NA               NA
## S_356         NA       NA        NA           NA         NA               NA
## S_438         NA       NA        NA           NA         NA               NA
## S_439         NA       NA        NA           NA         NA               NA
## S_440         NA       NA        NA           NA         NA               NA
## S_441         NA       NA        NA           NA         NA               NA
## S_442         NA       NA        NA           NA         NA               NA
## S_443         NA       NA        NA           NA         NA               NA
## S_444         NA       NA        NA           NA         NA               NA
## S_445         NA       NA        NA           NA         NA               NA
## S_446         NA       NA        NA           NA         NA               NA
## S_447         NA       NA        NA           NA         NA               NA
## S_448         NA       NA        NA           NA         NA               NA
## S_449         NA       NA        NA           NA         NA               NA
## S_450         NA       NA        NA           NA         NA               NA
## S_465         NA       NA        NA           NA         NA               NA
## S_466         NA       NA        NA           NA         NA               NA
## S_467         NA       NA        NA           NA         NA               NA
## S_468         NA       NA        NA           NA         NA               NA
## S_469         NA       NA        NA           NA         NA               NA
## S_470         NA       NA        NA           NA         NA               NA
## S_471         NA       NA        NA           NA         NA               NA
## S_472         NA       NA        NA           NA         NA               NA
## S_473         NA       NA        NA           NA         NA               NA
## S_474         NA       NA        NA           NA         NA               NA
## S_475         NA       NA        NA           NA         NA               NA
## S_476         NA       NA        NA           NA         NA               NA
## S_477         NA       NA        NA           NA         NA               NA
## S_478         NA       NA        NA           NA         NA               NA
##       QPCR_Enterobacteriaceae QPCR_Lachnospiraceae QPCR_Ruminococcaceae
## S_219                      NA                   NA                   NA
## S_220                      NA                   NA                   NA
## S_221                      NA                   NA                   NA
## S_222                      NA                   NA                   NA
## S_223                      NA                   NA                   NA
## S_224                      NA                   NA                   NA
## S_225                      NA                   NA                   NA
## S_226                      NA                   NA                   NA
## S_227                      NA                   NA                   NA
## S_228                      NA                   NA                   NA
## S_229                      NA                   NA                   NA
## S_230                      NA                   NA                   NA
## S_231                      NA                   NA                   NA
## S_232                      NA                   NA                   NA
## S_233                      NA                   NA                   NA
## S_234                      NA                   NA                   NA
## S_251                      NA                   NA                   NA
## S_252                      NA                   NA                   NA
## S_253                      NA                   NA                   NA
## S_254                      NA                   NA                   NA
## S_255                      NA                   NA                   NA
## S_256                      NA                   NA                   NA
## S_257                      NA                   NA                   NA
## S_259                      NA                   NA                   NA
## S_260                      NA                   NA                   NA
## S_261                      NA                   NA                   NA
## S_262                      NA                   NA                   NA
## S_263                      NA                   NA                   NA
## S_264                      NA                   NA                   NA
## S_265                      NA                   NA                   NA
## S_267                      NA                   NA                   NA
## S_268                      NA                   NA                   NA
## S_269                      NA                   NA                   NA
## S_270                      NA                   NA                   NA
## S_271                      NA                   NA                   NA
## S_272                      NA                   NA                   NA
## S_273                      NA                   NA                   NA
## S_274                      NA                   NA                   NA
## S_275                      NA                   NA                   NA
## S_276                      NA                   NA                   NA
## S_277                      NA                   NA                   NA
## S_278                      NA                   NA                   NA
## S_279                      NA                   NA                   NA
## S_280                      NA                   NA                   NA
## S_281                      NA                   NA                   NA
## S_298                      NA                   NA                   NA
## S_299                      NA                   NA                   NA
## S_300                      NA                   NA                   NA
## S_301                      NA                   NA                   NA
## S_302                      NA                   NA                   NA
## S_303                      NA                   NA                   NA
## S_304                      NA                   NA                   NA
## S_305                      NA                   NA                   NA
## S_306                      NA                   NA                   NA
## S_308                      NA                   NA                   NA
## S_309                      NA                   NA                   NA
## S_310                      NA                   NA                   NA
## S_311                      NA                   NA                   NA
## S_312                      NA                   NA                   NA
## S_313                      NA                   NA                   NA
## S_315                      NA                   NA                   NA
## S_316                      NA                   NA                   NA
## S_317                      NA                   NA                   NA
## S_318                      NA                   NA                   NA
## S_319                      NA                   NA                   NA
## S_320                      NA                   NA                   NA
## S_321                      NA                   NA                   NA
## S_322                      NA                   NA                   NA
## S_323                      NA                   NA                   NA
## S_324                      NA                   NA                   NA
## S_325                      NA                   NA                   NA
## S_326                      NA                   NA                   NA
## S_327                      NA                   NA                   NA
## S_341                      NA                   NA                   NA
## S_342                      NA                   NA                   NA
## S_343                      NA                   NA                   NA
## S_344                      NA                   NA                   NA
## S_345                      NA                   NA                   NA
## S_346                      NA                   NA                   NA
## S_347                      NA                   NA                   NA
## S_348                      NA                   NA                   NA
## S_349                      NA                   NA                   NA
## S_350                      NA                   NA                   NA
## S_351                      NA                   NA                   NA
## S_352                      NA                   NA                   NA
## S_353                      NA                   NA                   NA
## S_354                      NA                   NA                   NA
## S_355                      NA                   NA                   NA
## S_356                      NA                   NA                   NA
## S_438                      NA                   NA                   NA
## S_439                      NA                   NA                   NA
## S_440                      NA                   NA                   NA
## S_441                      NA                   NA                   NA
## S_442                      NA                   NA                   NA
## S_443                      NA                   NA                   NA
## S_444                      NA                   NA                   NA
## S_445                      NA                   NA                   NA
## S_446                      NA                   NA                   NA
## S_447                      NA                   NA                   NA
## S_448                      NA                   NA                   NA
## S_449                      NA                   NA                   NA
## S_450                      NA                   NA                   NA
## S_465                      NA                   NA                   NA
## S_466                      NA                   NA                   NA
## S_467                      NA                   NA                   NA
## S_468                      NA                   NA                   NA
## S_469                      NA                   NA                   NA
## S_470                      NA                   NA                   NA
## S_471                      NA                   NA                   NA
## S_472                      NA                   NA                   NA
## S_473                      NA                   NA                   NA
## S_474                      NA                   NA                   NA
## S_475                      NA                   NA                   NA
## S_476                      NA                   NA                   NA
## S_477                      NA                   NA                   NA
## S_478                      NA                   NA                   NA
##       QPCR_ButCoA Tyramine Putrescine Cadaverine Spermidine Spermine
## S_219          NA       NA         NA         NA         NA       NA
## S_220          NA       NA         NA         NA         NA       NA
## S_221          NA       NA         NA         NA         NA       NA
## S_222          NA       NA         NA         NA         NA       NA
## S_223          NA       NA         NA         NA         NA       NA
## S_224          NA       NA         NA         NA         NA       NA
## S_225          NA       NA         NA         NA         NA       NA
## S_226          NA       NA         NA         NA         NA       NA
## S_227          NA       NA         NA         NA         NA       NA
## S_228          NA       NA         NA         NA         NA       NA
## S_229          NA       NA         NA         NA         NA       NA
## S_230          NA       NA         NA         NA         NA       NA
## S_231          NA       NA         NA         NA         NA       NA
## S_232          NA       NA         NA         NA         NA       NA
## S_233          NA       NA         NA         NA         NA       NA
## S_234          NA       NA         NA         NA         NA       NA
## S_251          NA       NA         NA         NA         NA       NA
## S_252          NA       NA         NA         NA         NA       NA
## S_253          NA       NA         NA         NA         NA       NA
## S_254          NA       NA         NA         NA         NA       NA
## S_255          NA       NA         NA         NA         NA       NA
## S_256          NA       NA         NA         NA         NA       NA
## S_257          NA       NA         NA         NA         NA       NA
## S_259          NA       NA         NA         NA         NA       NA
## S_260          NA       NA         NA         NA         NA       NA
## S_261          NA       NA         NA         NA         NA       NA
## S_262          NA       NA         NA         NA         NA       NA
## S_263          NA       NA         NA         NA         NA       NA
## S_264          NA       NA         NA         NA         NA       NA
## S_265          NA       NA         NA         NA         NA       NA
## S_267          NA       NA         NA         NA         NA       NA
## S_268          NA       NA         NA         NA         NA       NA
## S_269          NA       NA         NA         NA         NA       NA
## S_270          NA       NA         NA         NA         NA       NA
## S_271          NA       NA         NA         NA         NA       NA
## S_272          NA       NA         NA         NA         NA       NA
## S_273          NA       NA         NA         NA         NA       NA
## S_274          NA       NA         NA         NA         NA       NA
## S_275          NA       NA         NA         NA         NA       NA
## S_276          NA       NA         NA         NA         NA       NA
## S_277          NA       NA         NA         NA         NA       NA
## S_278          NA       NA         NA         NA         NA       NA
## S_279          NA       NA         NA         NA         NA       NA
## S_280          NA       NA         NA         NA         NA       NA
## S_281          NA       NA         NA         NA         NA       NA
## S_298          NA       NA         NA         NA         NA       NA
## S_299          NA       NA         NA         NA         NA       NA
## S_300          NA       NA         NA         NA         NA       NA
## S_301          NA       NA         NA         NA         NA       NA
## S_302          NA       NA         NA         NA         NA       NA
## S_303          NA       NA         NA         NA         NA       NA
## S_304          NA       NA         NA         NA         NA       NA
## S_305          NA       NA         NA         NA         NA       NA
## S_306          NA       NA         NA         NA         NA       NA
## S_308          NA       NA         NA         NA         NA       NA
## S_309          NA       NA         NA         NA         NA       NA
## S_310          NA       NA         NA         NA         NA       NA
## S_311          NA       NA         NA         NA         NA       NA
## S_312          NA       NA         NA         NA         NA       NA
## S_313          NA       NA         NA         NA         NA       NA
## S_315          NA       NA         NA         NA         NA       NA
## S_316          NA       NA         NA         NA         NA       NA
## S_317          NA       NA         NA         NA         NA       NA
## S_318          NA       NA         NA         NA         NA       NA
## S_319          NA       NA         NA         NA         NA       NA
## S_320          NA       NA         NA         NA         NA       NA
## S_321          NA       NA         NA         NA         NA       NA
## S_322          NA       NA         NA         NA         NA       NA
## S_323          NA       NA         NA         NA         NA       NA
## S_324          NA       NA         NA         NA         NA       NA
## S_325          NA       NA         NA         NA         NA       NA
## S_326          NA       NA         NA         NA         NA       NA
## S_327          NA       NA         NA         NA         NA       NA
## S_341          NA       NA         NA         NA         NA       NA
## S_342          NA       NA         NA         NA         NA       NA
## S_343          NA       NA         NA         NA         NA       NA
## S_344          NA       NA         NA         NA         NA       NA
## S_345          NA       NA         NA         NA         NA       NA
## S_346          NA       NA         NA         NA         NA       NA
## S_347          NA       NA         NA         NA         NA       NA
## S_348          NA       NA         NA         NA         NA       NA
## S_349          NA       NA         NA         NA         NA       NA
## S_350          NA       NA         NA         NA         NA       NA
## S_351          NA       NA         NA         NA         NA       NA
## S_352          NA       NA         NA         NA         NA       NA
## S_353          NA       NA         NA         NA         NA       NA
## S_354          NA       NA         NA         NA         NA       NA
## S_355          NA       NA         NA         NA         NA       NA
## S_356          NA       NA         NA         NA         NA       NA
## S_438          NA       NA         NA         NA         NA       NA
## S_439          NA       NA         NA         NA         NA       NA
## S_440          NA       NA         NA         NA         NA       NA
## S_441          NA       NA         NA         NA         NA       NA
## S_442          NA       NA         NA         NA         NA       NA
## S_443          NA       NA         NA         NA         NA       NA
## S_444          NA       NA         NA         NA         NA       NA
## S_445          NA       NA         NA         NA         NA       NA
## S_446          NA       NA         NA         NA         NA       NA
## S_447          NA       NA         NA         NA         NA       NA
## S_448          NA       NA         NA         NA         NA       NA
## S_449          NA       NA         NA         NA         NA       NA
## S_450          NA       NA         NA         NA         NA       NA
## S_465          NA       NA         NA         NA         NA       NA
## S_466          NA       NA         NA         NA         NA       NA
## S_467          NA       NA         NA         NA         NA       NA
## S_468          NA       NA         NA         NA         NA       NA
## S_469          NA       NA         NA         NA         NA       NA
## S_470          NA       NA         NA         NA         NA       NA
## S_471          NA       NA         NA         NA         NA       NA
## S_472          NA       NA         NA         NA         NA       NA
## S_473          NA       NA         NA         NA         NA       NA
## S_474          NA       NA         NA         NA         NA       NA
## S_475          NA       NA         NA         NA         NA       NA
## S_476          NA       NA         NA         NA         NA       NA
## S_477          NA       NA         NA         NA         NA       NA
## S_478          NA       NA         NA         NA         NA       NA
##       B24_Acetate B24_Propionate B24_Butyrate B24_Formate B24_Succinate
## S_219          NA             NA           NA          NA            NA
## S_220          NA             NA           NA          NA            NA
## S_221          NA             NA           NA          NA            NA
## S_222          NA             NA           NA          NA            NA
## S_223          NA             NA           NA          NA            NA
## S_224          NA             NA           NA          NA            NA
## S_225          NA             NA           NA          NA            NA
## S_226          NA             NA           NA          NA            NA
## S_227          NA             NA           NA          NA            NA
## S_228          NA             NA           NA          NA            NA
## S_229          NA             NA           NA          NA            NA
## S_230          NA             NA           NA          NA            NA
## S_231          NA             NA           NA          NA            NA
## S_232          NA             NA           NA          NA            NA
## S_233          NA             NA           NA          NA            NA
## S_234          NA             NA           NA          NA            NA
## S_251          NA             NA           NA          NA            NA
## S_252          NA             NA           NA          NA            NA
## S_253          NA             NA           NA          NA            NA
## S_254          NA             NA           NA          NA            NA
## S_255          NA             NA           NA          NA            NA
## S_256          NA             NA           NA          NA            NA
## S_257          NA             NA           NA          NA            NA
## S_259          NA             NA           NA          NA            NA
## S_260          NA             NA           NA          NA            NA
## S_261          NA             NA           NA          NA            NA
## S_262          NA             NA           NA          NA            NA
## S_263          NA             NA           NA          NA            NA
## S_264          NA             NA           NA          NA            NA
## S_265          NA             NA           NA          NA            NA
## S_267          NA             NA           NA          NA            NA
## S_268          NA             NA           NA          NA            NA
## S_269          NA             NA           NA          NA            NA
## S_270          NA             NA           NA          NA            NA
## S_271          NA             NA           NA          NA            NA
## S_272          NA             NA           NA          NA            NA
## S_273          NA             NA           NA          NA            NA
## S_274          NA             NA           NA          NA            NA
## S_275          NA             NA           NA          NA            NA
## S_276          NA             NA           NA          NA            NA
## S_277          NA             NA           NA          NA            NA
## S_278          NA             NA           NA          NA            NA
## S_279          NA             NA           NA          NA            NA
## S_280          NA             NA           NA          NA            NA
## S_281          NA             NA           NA          NA            NA
## S_298          NA             NA           NA          NA            NA
## S_299          NA             NA           NA          NA            NA
## S_300          NA             NA           NA          NA            NA
## S_301          NA             NA           NA          NA            NA
## S_302          NA             NA           NA          NA            NA
## S_303          NA             NA           NA          NA            NA
## S_304          NA             NA           NA          NA            NA
## S_305          NA             NA           NA          NA            NA
## S_306          NA             NA           NA          NA            NA
## S_308          NA             NA           NA          NA            NA
## S_309          NA             NA           NA          NA            NA
## S_310          NA             NA           NA          NA            NA
## S_311          NA             NA           NA          NA            NA
## S_312          NA             NA           NA          NA            NA
## S_313          NA             NA           NA          NA            NA
## S_315          NA             NA           NA          NA            NA
## S_316          NA             NA           NA          NA            NA
## S_317          NA             NA           NA          NA            NA
## S_318          NA             NA           NA          NA            NA
## S_319          NA             NA           NA          NA            NA
## S_320          NA             NA           NA          NA            NA
## S_321          NA             NA           NA          NA            NA
## S_322          NA             NA           NA          NA            NA
## S_323          NA             NA           NA          NA            NA
## S_324          NA             NA           NA          NA            NA
## S_325          NA             NA           NA          NA            NA
## S_326          NA             NA           NA          NA            NA
## S_327          NA             NA           NA          NA            NA
## S_341          NA             NA           NA          NA            NA
## S_342          NA             NA           NA          NA            NA
## S_343          NA             NA           NA          NA            NA
## S_344          NA             NA           NA          NA            NA
## S_345          NA             NA           NA          NA            NA
## S_346          NA             NA           NA          NA            NA
## S_347          NA             NA           NA          NA            NA
## S_348          NA             NA           NA          NA            NA
## S_349          NA             NA           NA          NA            NA
## S_350          NA             NA           NA          NA            NA
## S_351          NA             NA           NA          NA            NA
## S_352          NA             NA           NA          NA            NA
## S_353          NA             NA           NA          NA            NA
## S_354          NA             NA           NA          NA            NA
## S_355          NA             NA           NA          NA            NA
## S_356          NA             NA           NA          NA            NA
## S_438          NA             NA           NA          NA            NA
## S_439          NA             NA           NA          NA            NA
## S_440          NA             NA           NA          NA            NA
## S_441          NA             NA           NA          NA            NA
## S_442          NA             NA           NA          NA            NA
## S_443          NA             NA           NA          NA            NA
## S_444          NA             NA           NA          NA            NA
## S_445          NA             NA           NA          NA            NA
## S_446          NA             NA           NA          NA            NA
## S_447          NA             NA           NA          NA            NA
## S_448          NA             NA           NA          NA            NA
## S_449          NA             NA           NA          NA            NA
## S_450          NA             NA           NA          NA            NA
## S_465          NA             NA           NA          NA            NA
## S_466          NA             NA           NA          NA            NA
## S_467          NA             NA           NA          NA            NA
## S_468          NA             NA           NA          NA            NA
## S_469          NA             NA           NA          NA            NA
## S_470          NA             NA           NA          NA            NA
## S_471          NA             NA           NA          NA            NA
## S_472          NA             NA           NA          NA            NA
## S_473          NA             NA           NA          NA            NA
## S_474          NA             NA           NA          NA            NA
## S_475          NA             NA           NA          NA            NA
## S_476          NA             NA           NA          NA            NA
## S_477          NA             NA           NA          NA            NA
## S_478          NA             NA           NA          NA            NA
##       B24_Lactate B24_BCFA succinate lactate formate acetate propionate
## S_219          NA       NA        NA      NA      NA      NA         NA
## S_220          NA       NA        NA      NA      NA      NA         NA
## S_221          NA       NA        NA      NA      NA      NA         NA
## S_222          NA       NA        NA      NA      NA      NA         NA
## S_223          NA       NA        NA      NA      NA      NA         NA
## S_224          NA       NA        NA      NA      NA      NA         NA
## S_225          NA       NA        NA      NA      NA      NA         NA
## S_226          NA       NA        NA      NA      NA      NA         NA
## S_227          NA       NA        NA      NA      NA      NA         NA
## S_228          NA       NA        NA      NA      NA      NA         NA
## S_229          NA       NA        NA      NA      NA      NA         NA
## S_230          NA       NA        NA      NA      NA      NA         NA
## S_231          NA       NA        NA      NA      NA      NA         NA
## S_232          NA       NA        NA      NA      NA      NA         NA
## S_233          NA       NA        NA      NA      NA      NA         NA
## S_234          NA       NA        NA      NA      NA      NA         NA
## S_251          NA       NA        NA      NA      NA      NA         NA
## S_252          NA       NA        NA      NA      NA      NA         NA
## S_253          NA       NA        NA      NA      NA      NA         NA
## S_254          NA       NA        NA      NA      NA      NA         NA
## S_255          NA       NA        NA      NA      NA      NA         NA
## S_256          NA       NA        NA      NA      NA      NA         NA
## S_257          NA       NA        NA      NA      NA      NA         NA
## S_259          NA       NA        NA      NA      NA      NA         NA
## S_260          NA       NA        NA      NA      NA      NA         NA
## S_261          NA       NA        NA      NA      NA      NA         NA
## S_262          NA       NA        NA      NA      NA      NA         NA
## S_263          NA       NA        NA      NA      NA      NA         NA
## S_264          NA       NA        NA      NA      NA      NA         NA
## S_265          NA       NA        NA      NA      NA      NA         NA
## S_267          NA       NA        NA      NA      NA      NA         NA
## S_268          NA       NA        NA      NA      NA      NA         NA
## S_269          NA       NA        NA      NA      NA      NA         NA
## S_270          NA       NA        NA      NA      NA      NA         NA
## S_271          NA       NA        NA      NA      NA      NA         NA
## S_272          NA       NA        NA      NA      NA      NA         NA
## S_273          NA       NA        NA      NA      NA      NA         NA
## S_274          NA       NA        NA      NA      NA      NA         NA
## S_275          NA       NA        NA      NA      NA      NA         NA
## S_276          NA       NA        NA      NA      NA      NA         NA
## S_277          NA       NA        NA      NA      NA      NA         NA
## S_278          NA       NA        NA      NA      NA      NA         NA
## S_279          NA       NA        NA      NA      NA      NA         NA
## S_280          NA       NA        NA      NA      NA      NA         NA
## S_281          NA       NA        NA      NA      NA      NA         NA
## S_298          NA       NA        NA      NA      NA      NA         NA
## S_299          NA       NA        NA      NA      NA      NA         NA
## S_300          NA       NA        NA      NA      NA      NA         NA
## S_301          NA       NA        NA      NA      NA      NA         NA
## S_302          NA       NA        NA      NA      NA      NA         NA
## S_303          NA       NA        NA      NA      NA      NA         NA
## S_304          NA       NA        NA      NA      NA      NA         NA
## S_305          NA       NA        NA      NA      NA      NA         NA
## S_306          NA       NA        NA      NA      NA      NA         NA
## S_308          NA       NA        NA      NA      NA      NA         NA
## S_309          NA       NA        NA      NA      NA      NA         NA
## S_310          NA       NA        NA      NA      NA      NA         NA
## S_311          NA       NA        NA      NA      NA      NA         NA
## S_312          NA       NA        NA      NA      NA      NA         NA
## S_313          NA       NA        NA      NA      NA      NA         NA
## S_315          NA       NA        NA      NA      NA      NA         NA
## S_316          NA       NA        NA      NA      NA      NA         NA
## S_317          NA       NA        NA      NA      NA      NA         NA
## S_318          NA       NA        NA      NA      NA      NA         NA
## S_319          NA       NA        NA      NA      NA      NA         NA
## S_320          NA       NA        NA      NA      NA      NA         NA
## S_321          NA       NA        NA      NA      NA      NA         NA
## S_322          NA       NA        NA      NA      NA      NA         NA
## S_323          NA       NA        NA      NA      NA      NA         NA
## S_324          NA       NA        NA      NA      NA      NA         NA
## S_325          NA       NA        NA      NA      NA      NA         NA
## S_326          NA       NA        NA      NA      NA      NA         NA
## S_327          NA       NA        NA      NA      NA      NA         NA
## S_341          NA       NA        NA      NA      NA      NA         NA
## S_342          NA       NA        NA      NA      NA      NA         NA
## S_343          NA       NA        NA      NA      NA      NA         NA
## S_344          NA       NA        NA      NA      NA      NA         NA
## S_345          NA       NA        NA      NA      NA      NA         NA
## S_346          NA       NA        NA      NA      NA      NA         NA
## S_347          NA       NA        NA      NA      NA      NA         NA
## S_348          NA       NA        NA      NA      NA      NA         NA
## S_349          NA       NA        NA      NA      NA      NA         NA
## S_350          NA       NA        NA      NA      NA      NA         NA
## S_351          NA       NA        NA      NA      NA      NA         NA
## S_352          NA       NA        NA      NA      NA      NA         NA
## S_353          NA       NA        NA      NA      NA      NA         NA
## S_354          NA       NA        NA      NA      NA      NA         NA
## S_355          NA       NA        NA      NA      NA      NA         NA
## S_356          NA       NA        NA      NA      NA      NA         NA
## S_438          NA       NA        NA      NA      NA      NA         NA
## S_439          NA       NA        NA      NA      NA      NA         NA
## S_440          NA       NA        NA      NA      NA      NA         NA
## S_441          NA       NA        NA      NA      NA      NA         NA
## S_442          NA       NA        NA      NA      NA      NA         NA
## S_443          NA       NA        NA      NA      NA      NA         NA
## S_444          NA       NA        NA      NA      NA      NA         NA
## S_445          NA       NA        NA      NA      NA      NA         NA
## S_446          NA       NA        NA      NA      NA      NA         NA
## S_447          NA       NA        NA      NA      NA      NA         NA
## S_448          NA       NA        NA      NA      NA      NA         NA
## S_449          NA       NA        NA      NA      NA      NA         NA
## S_450          NA       NA        NA      NA      NA      NA         NA
## S_465          NA       NA        NA      NA      NA      NA         NA
## S_466          NA       NA        NA      NA      NA      NA         NA
## S_467          NA       NA        NA      NA      NA      NA         NA
## S_468          NA       NA        NA      NA      NA      NA         NA
## S_469          NA       NA        NA      NA      NA      NA         NA
## S_470          NA       NA        NA      NA      NA      NA         NA
## S_471          NA       NA        NA      NA      NA      NA         NA
## S_472          NA       NA        NA      NA      NA      NA         NA
## S_473          NA       NA        NA      NA      NA      NA         NA
## S_474          NA       NA        NA      NA      NA      NA         NA
## S_475          NA       NA        NA      NA      NA      NA         NA
## S_476          NA       NA        NA      NA      NA      NA         NA
## S_477          NA       NA        NA      NA      NA      NA         NA
## S_478          NA       NA        NA      NA      NA      NA         NA
##       isobutyrate butyrate isovalerate valerate BCFA total is.neg period
## S_219          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_220          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_221          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_222          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_223          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_224          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_225          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_226          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_227          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_228          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_229          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_230          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_231          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_232          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_233          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_234          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_251          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_252          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_253          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_254          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_255          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_256          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_257          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_259          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_260          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_261          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_262          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_263          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_264          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_265          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_267          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_268          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_269          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_270          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_271          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_272          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_273          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_274          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_275          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_276          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_277          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_278          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_279          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_280          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_281          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_298          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_299          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_300          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_301          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_302          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_303          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_304          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_305          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_306          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_308          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_309          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_310          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_311          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_312          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_313          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_315          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_316          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_317          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_318          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_319          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_320          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_321          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_322          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_323          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_324          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_325          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_326          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_327          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_341          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_342          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_343          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_344          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_345          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_346          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_347          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_348          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_349          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_350          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_351          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_352          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_353          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_354          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_355          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_356          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_438          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_439          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_440          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_441          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_442          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_443          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_444          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_445          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_446          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_447          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_448          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_449          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_450          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_465          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_466          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_467          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_468          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_469          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_470          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_471          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_472          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_473          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_474          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_475          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_476          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_477          NA       NA          NA       NA   NA    NA  FALSE    pNA
## S_478          NA       NA          NA       NA   NA    NA  FALSE    pNA
##       day_fermentation_d IR_TR Stab_Treat experiment_period period_treat
## S_219                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_220                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_221                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_222                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_223                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_224                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_225                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_226                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_227                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_228                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_229                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_230                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_231                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_232                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_233                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_234                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_251                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_252                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_253                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_254                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_255                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_256                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_257                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_259                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_260                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_261                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_262                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_263                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_264                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_265                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_267                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_268                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_269                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_270                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_271                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_272                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_273                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_274                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_275                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_276                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_277                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_278                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_279                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_280                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_281                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_298                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_299                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_300                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_301                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_302                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_303                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_304                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_305                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_306                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_308                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_309                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_310                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_311                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_312                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_313                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_315                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_316                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_317                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_318                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_319                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_320                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_321                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_322                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_323                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_324                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_325                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_326                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_327                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_341                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_342                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_343                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_344                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_345                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_346                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_347                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_348                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_349                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_350                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_351                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_352                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_353                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_354                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_355                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_356                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_438                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_439                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_440                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_441                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_442                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_443                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_444                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_445                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_446                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_447                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_448                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_449                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_450                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_465                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_466                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_467                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_468                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_469                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_470                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_471                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_472                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_473                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_474                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_475                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_476                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_477                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
## S_478                 NA  <NA>       <NA>         mIMT_2_NA         <NA>
##       period_reactor treatment_grouped treatment_invivo day_num
## S_219         pNA_NA          Eubiotic             none      -1
## S_220         pNA_NA          Eubiotic             none      -1
## S_221         pNA_NA          Eubiotic             none      -1
## S_222         pNA_NA          Eubiotic             none      -1
## S_223         pNA_NA          Eubiotic             none      -1
## S_224         pNA_NA          Eubiotic             none      -1
## S_225         pNA_NA          Eubiotic             none      -1
## S_226         pNA_NA          Eubiotic             none      -1
## S_227         pNA_NA         Dysbiotic             none      -1
## S_228         pNA_NA         Dysbiotic             none      -1
## S_229         pNA_NA         Dysbiotic             none      -1
## S_230         pNA_NA         Dysbiotic             none      -1
## S_231         pNA_NA         Dysbiotic             none      -1
## S_232         pNA_NA         Dysbiotic             none      -1
## S_233         pNA_NA         Dysbiotic             none      -1
## S_234         pNA_NA         Dysbiotic             none      -1
## S_251         pNA_NA               DSS             none      -1
## S_252         pNA_NA               DSS             none      -1
## S_253         pNA_NA               DSS             none      -1
## S_254         pNA_NA               DSS             none      -1
## S_255         pNA_NA               DSS             none      -1
## S_256         pNA_NA               DSS             none      -1
## S_257         pNA_NA               DSS             none      -1
## S_259         pNA_NA               H2O             none      -1
## S_260         pNA_NA               H2O             none      -1
## S_261         pNA_NA               H2O             none      -1
## S_262         pNA_NA               H2O             none      -1
## S_263         pNA_NA               H2O             none      -1
## S_264         pNA_NA               H2O             none      -1
## S_265         pNA_NA               H2O             none      -1
## S_267         pNA_NA          Eubiotic         Eubiotic      12
## S_268         pNA_NA          Eubiotic         Eubiotic      12
## S_269         pNA_NA          Eubiotic         Eubiotic      12
## S_270         pNA_NA          Eubiotic         Eubiotic      12
## S_271         pNA_NA          Eubiotic         Eubiotic      12
## S_272         pNA_NA          Eubiotic         Eubiotic      12
## S_273         pNA_NA          Eubiotic         Eubiotic      12
## S_274         pNA_NA         Dysbiotic        Dysbiotic      12
## S_275         pNA_NA         Dysbiotic        Dysbiotic      12
## S_276         pNA_NA         Dysbiotic        Dysbiotic      12
## S_277         pNA_NA         Dysbiotic        Dysbiotic      12
## S_278         pNA_NA         Dysbiotic        Dysbiotic      12
## S_279         pNA_NA         Dysbiotic        Dysbiotic      12
## S_280         pNA_NA         Dysbiotic        Dysbiotic      12
## S_281         pNA_NA         Dysbiotic        Dysbiotic      12
## S_298         pNA_NA               DSS              DSS      12
## S_299         pNA_NA               DSS              DSS      12
## S_300         pNA_NA               DSS              DSS      12
## S_301         pNA_NA               DSS              DSS      12
## S_302         pNA_NA               DSS              DSS      12
## S_303         pNA_NA               DSS              DSS      12
## S_304         pNA_NA               DSS              DSS      12
## S_305         pNA_NA               H2O             none      12
## S_306         pNA_NA               H2O             none      12
## S_308         pNA_NA               H2O             none      12
## S_309         pNA_NA               H2O             none      12
## S_310         pNA_NA               H2O             none      12
## S_311         pNA_NA               H2O             none      12
## S_312         pNA_NA               H2O             none      12
## S_313         pNA_NA          Eubiotic         Eubiotic      14
## S_315         pNA_NA          Eubiotic         Eubiotic      14
## S_316         pNA_NA          Eubiotic         Eubiotic      14
## S_317         pNA_NA          Eubiotic         Eubiotic      14
## S_318         pNA_NA          Eubiotic         Eubiotic      14
## S_319         pNA_NA          Eubiotic         Eubiotic      14
## S_320         pNA_NA         Dysbiotic        Dysbiotic      14
## S_321         pNA_NA         Dysbiotic        Dysbiotic      14
## S_322         pNA_NA         Dysbiotic        Dysbiotic      14
## S_323         pNA_NA         Dysbiotic        Dysbiotic      14
## S_324         pNA_NA         Dysbiotic        Dysbiotic      14
## S_325         pNA_NA         Dysbiotic        Dysbiotic      14
## S_326         pNA_NA         Dysbiotic        Dysbiotic      14
## S_327         pNA_NA         Dysbiotic        Dysbiotic      14
## S_341         pNA_NA               DSS              DSS      14
## S_342         pNA_NA               DSS              DSS      14
## S_343         pNA_NA               DSS              DSS      14
## S_344         pNA_NA               DSS              DSS      14
## S_345         pNA_NA               DSS              DSS      14
## S_346         pNA_NA               DSS              DSS      14
## S_347         pNA_NA               DSS              DSS      14
## S_348         pNA_NA               DSS              DSS      14
## S_349         pNA_NA               H2O             none      14
## S_350         pNA_NA               H2O             none      14
## S_351         pNA_NA               H2O             none      14
## S_352         pNA_NA               H2O             none      14
## S_353         pNA_NA               H2O             none      14
## S_354         pNA_NA               H2O             none      14
## S_355         pNA_NA               H2O             none      14
## S_356         pNA_NA               H2O             none      14
## S_438         pNA_NA          Eubiotic         Eubiotic       6
## S_439         pNA_NA          Eubiotic         Eubiotic       6
## S_440         pNA_NA          Eubiotic         Eubiotic       6
## S_441         pNA_NA          Eubiotic         Eubiotic       6
## S_442         pNA_NA          Eubiotic         Eubiotic       6
## S_443         pNA_NA         Dysbiotic        Dysbiotic       6
## S_444         pNA_NA         Dysbiotic        Dysbiotic       6
## S_445         pNA_NA         Dysbiotic        Dysbiotic       6
## S_446         pNA_NA         Dysbiotic        Dysbiotic       6
## S_447         pNA_NA         Dysbiotic        Dysbiotic       6
## S_448         pNA_NA         Dysbiotic        Dysbiotic       6
## S_449         pNA_NA         Dysbiotic        Dysbiotic       6
## S_450         pNA_NA         Dysbiotic        Dysbiotic       6
## S_465         pNA_NA               DSS              DSS       6
## S_466         pNA_NA               DSS              DSS       6
## S_467         pNA_NA               DSS              DSS       6
## S_468         pNA_NA               DSS              DSS       6
## S_469         pNA_NA               DSS              DSS       6
## S_470         pNA_NA               DSS              DSS       6
## S_471         pNA_NA               H2O             none       6
## S_472         pNA_NA               H2O             none       6
## S_473         pNA_NA               H2O             none       6
## S_474         pNA_NA               H2O             none       6
## S_475         pNA_NA               H2O             none       6
## S_476         pNA_NA               H2O             none       6
## S_477         pNA_NA               H2O             none       6
## S_478         pNA_NA               H2O             none       6
```


```r
ps_invivo %>% 
  sample_data() %>% 
  data.frame() %>% 
  select(-reads ,-remove,-experiment, -owner, -run, -region_1, -responsible_seq,-facility, -sample_type, -reactor, -day_IR, -day_fermentation, -treatment_period, -treatment_exvivo, -experiment_period, -period_treat, -IR_TR, -period_treat, -sample_file, -colon_lenght:-period_reactor, -sample_label, -Fermentation, ) %>% 
  select(Sample, everything()) %>% 
  mutate(treatment_invivo = recode(treatment_invivo, 
                                   `H2O` = "Control",
                                   `Eubiotic` = "TreatA",
                                   `Dysbiotic` = "TreatB",
                                   `DSS` = "TreatC"),
         treatment = recode(treatment, 
                            `H2O` = "Control",
                            `Eubiotic` = "TreatA",
                            `Dysbiotic` = "TreatB",
                            `DSS` = "TreatC"),
         treatment_grouped = recode(treatment_grouped, 
                                    `H2O` = "Control" ,                
                                    `Eubiotic` = "TreatA",
                                    `Dysbiotic` = "TreatB",
                                    `DSS` = "TreatC")) -> sample_data(ps_invivo)
ps_invivo %>% 
  sample_data() %>% 
  data.frame()
```

```
##       Sample  input filtered filtered_pc denoisedF denoisedR denoisedF_pc
## S_219  S_219  44852    44233        0.99     42985     43681         0.97
## S_220  S_220  72545    72048        0.99     70726     70818         0.98
## S_221  S_221  76560    76051        0.99     74773     75012         0.98
## S_222  S_222  62479    62192        1.00     61528     61201         0.99
## S_223  S_223  78944    78670        1.00     77717     77643         0.99
## S_224  S_224  72930    72622        1.00     71294     71600         0.98
## S_225  S_225 194545   193572        0.99    191373    191090         0.99
## S_226  S_226  62034    61405        0.99     60035     60599         0.98
## S_227  S_227  77255    76602        0.99     75092     75664         0.98
## S_228  S_228  58784    58444        0.99     57469     57504         0.98
## S_229  S_229  65923    65701        1.00     64821     65102         0.99
## S_230  S_230  51065    50870        1.00     50385     50108         0.99
## S_231  S_231 123400   123134        1.00    122255    122251         0.99
## S_232  S_232  60612    60294        0.99     59542     59544         0.99
## S_233  S_233  66974    66745        1.00     65951     65893         0.99
## S_234  S_234  59225    58828        0.99     57814     57973         0.98
## S_251  S_251  46134    45771        0.99     44977     45491         0.98
## S_252  S_252  43200    42935        0.99     41926     42309         0.98
## S_253  S_253  61771    61392        0.99     60360     60789         0.98
## S_254  S_254  52064    51842        1.00     51293     51367         0.99
## S_255  S_255  60985    60720        1.00     60048     60142         0.99
## S_256  S_256  58725    58433        1.00     57619     57826         0.99
## S_257  S_257  53218    52990        1.00     52344     52408         0.99
## S_259  S_259  35191    34713        0.99     33684     34288         0.97
## S_260  S_260  59786    59402        0.99     58312     58638         0.98
## S_261  S_261  71389    70960        0.99     69815     70122         0.98
## S_262  S_262  48071    47899        1.00     47393     47356         0.99
## S_263  S_263  60459    60267        1.00     59414     59486         0.99
## S_264  S_264  56699    56368        0.99     55608     55822         0.99
## S_265  S_265  52415    52164        1.00     51488     51508         0.99
## S_267  S_267  40629    40017        0.98     37962     39383         0.95
## S_268  S_268  67522    67104        0.99     65052     65909         0.97
## S_269  S_269  62759    62323        0.99     60618     61516         0.97
## S_270  S_270  51923    51673        1.00     50793     51026         0.98
## S_271  S_271  67574    67308        1.00     66518     66567         0.99
## S_272  S_272  69843    69434        0.99     68180     68564         0.98
## S_273  S_273  53162    52908        1.00     52171     52277         0.99
## S_274  S_274  89399    88860        0.99     86777     87744         0.98
## S_275  S_275  47024    46664        0.99     46031     46150         0.99
## S_276  S_276  52772    52491        0.99     51730     51714         0.99
## S_277  S_277  57920    57707        1.00     56978     56801         0.99
## S_278  S_278  45769    45629        1.00     45263     44921         0.99
## S_279  S_279  69581    69393        1.00     68667     68504         0.99
## S_280  S_280  51654    51454        1.00     50897     50887         0.99
## S_281  S_281  47303    47147        1.00     46707     46496         0.99
## S_298  S_298  45576    45320        0.99     44563     44714         0.98
## S_299  S_299  37551    37094        0.99     36164     36676         0.97
## S_300  S_300  49047    48750        0.99     47852     48161         0.98
## S_301  S_301  55629    55307        0.99     54626     54653         0.99
## S_302  S_302  39112    38954        1.00     38658     38534         0.99
## S_303  S_303  56158    55984        1.00     54978     55559         0.98
## S_304  S_304  56174    55990        1.00     54558     55582         0.97
## S_305  S_305  49080    48860        1.00     48261     48353         0.99
## S_306  S_306  51145    50825        0.99     50082     50314         0.99
## S_308  S_308  61065    60799        1.00     59746     59875         0.98
## S_309  S_309  49604    49210        0.99     48411     48539         0.98
## S_310  S_310  39676    39568        1.00     39079     39070         0.99
## S_311  S_311  67509    67341        1.00     66594     66566         0.99
## S_312  S_312  54562    54327        1.00     53556     53619         0.99
## S_313  S_313  51990    51826        1.00     51343     51218         0.99
## S_315  S_315  56389    56103        0.99     55311     55282         0.99
## S_316  S_316  52602    52284        0.99     51592     51797         0.99
## S_317  S_317  55353    55045        0.99     54304     54490         0.99
## S_318  S_318  47228    47005        1.00     46369     46415         0.99
## S_319  S_319  51051    50825        1.00     50302     50300         0.99
## S_320  S_320  50339    50129        1.00     49527     49544         0.99
## S_321  S_321  53955    53733        1.00     53035     53076         0.99
## S_322  S_322  50561    50357        1.00     49892     49982         0.99
## S_323  S_323  55667    55424        1.00     54749     54623         0.99
## S_324  S_324  53070    52783        0.99     52078     52396         0.99
## S_325  S_325  43764    43485        0.99     42777     43048         0.98
## S_326  S_326  55187    55049        1.00     54210     54461         0.98
## S_327  S_327  44766    44572        1.00     44063     44115         0.99
## S_341  S_341  51618    51384        1.00     50750     50912         0.99
## S_342  S_342  38695    38576        1.00     38190     38209         0.99
## S_343  S_343  31639    31476        0.99     31038     31086         0.99
## S_344  S_344  57314    57167        1.00     56524     56611         0.99
## S_345  S_345  26796    26719        1.00     26323     26354         0.99
## S_346  S_346  42867    42762        1.00     42236     42312         0.99
## S_347  S_347  49591    49374        1.00     48708     48943         0.99
## S_348  S_348  49706    49595        1.00     49100     49181         0.99
## S_349  S_349  82068    81706        1.00     80665     80976         0.99
## S_350  S_350  52750    52500        1.00     51830     51938         0.99
## S_351  S_351  73446    73177        1.00     72467     72411         0.99
## S_352  S_352  69764    69491        1.00     68538     68691         0.99
## S_353  S_353  61764    61575        1.00     60824     61048         0.99
## S_354  S_354  52313    52053        1.00     51372     51569         0.99
## S_355  S_355  51593    51313        0.99     50391     50675         0.98
## S_356  S_356  46806    46542        0.99     45968     45972         0.99
## S_438  S_438  66744    66455        1.00     65323     65479         0.98
## S_439  S_439  65688    65409        1.00     64553     64591         0.99
## S_440  S_440  42282    42132        1.00     41655     41507         0.99
## S_441  S_441  70679    70422        1.00     69786     69632         0.99
## S_442  S_442  64149    63878        1.00     62852     63031         0.98
## S_443  S_443  66106    65772        0.99     65139     65014         0.99
## S_444  S_444  67344    67025        1.00     65804     66358         0.98
## S_445  S_445  47094    46707        0.99     45473     46244         0.97
## S_446  S_446  67572    67267        1.00     65981     66428         0.98
## S_447  S_447  67255    66994        1.00     65937     66416         0.98
## S_448  S_448  57626    57445        1.00     56882     56862         0.99
## S_449  S_449  69800    69629        1.00     69042     68944         0.99
## S_450  S_450  66731    66475        1.00     65494     65868         0.99
## S_465  S_465  69793    69596        1.00     68933     69150         0.99
## S_466  S_466  73792    73415        0.99     72195     72677         0.98
## S_467  S_467  64707    64401        1.00     63506     63635         0.99
## S_468  S_468  69869    69411        0.99     67836     68545         0.98
## S_469  S_469  22876    22575        0.99     21552     22319         0.95
## S_470  S_470  70853    70436        0.99     69010     69691         0.98
## S_471  S_471  66286    66030        1.00     64938     65385         0.98
## S_472  S_472  63794    63539        1.00     62844     62818         0.99
## S_473  S_473  76951    76627        1.00     76016     76028         0.99
## S_474  S_474  61480    61249        1.00     60023     60475         0.98
## S_475  S_475  65287    65039        1.00     64261     64598         0.99
## S_476  S_476  62397    61982        0.99     60877     61304         0.98
## S_477  S_477  40902    40453        0.99     39440     40103         0.97
## S_478  S_478  67854    67566        1.00     66284     66723         0.98
##       denoisedR_pc merged merged_pc tabled chimera_out length_filtered
## S_219         0.99  37694      0.88  37694       26929           26929
## S_220         0.98  63881      0.90  63881       44125           44125
## S_221         0.99  67904      0.91  67904       49814           49814
## S_222         0.98  54905      0.89  54905       42532           42532
## S_223         0.99  68254      0.88  68254       49859           49859
## S_224         0.99  60502      0.85  60502       44187           44187
## S_225         0.99 173270      0.91 173270      124781          124781
## S_226         0.99  53780      0.90  53780       39110           39110
## S_227         0.99  65965      0.88  65965       47209           47209
## S_228         0.98  52534      0.91  52534       39238           39238
## S_229         0.99  59729      0.92  59729       49569           49569
## S_230         0.99  45226      0.90  45226       32999           32999
## S_231         0.99 113946      0.93 113946       89836           89836
## S_232         0.99  53809      0.90  53809       37806           37806
## S_233         0.99  58534      0.89  58534       44726           44726
## S_234         0.99  52635      0.91  52635       37761           37761
## S_251         0.99  42248      0.94  42248       36121           36121
## S_252         0.99  37345      0.89  37345       27664           27664
## S_253         0.99  55612      0.92  55612       39123           39123
## S_254         0.99  47517      0.93  47517       34414           34414
## S_255         0.99  55956      0.93  55956       41213           41213
## S_256         0.99  53508      0.93  53508       39884           39884
## S_257         0.99  48415      0.92  48415       35603           35603
## S_259         0.99  29649      0.88  29649       22393           22393
## S_260         0.99  53165      0.91  53165       37401           37401
## S_261         0.99  63752      0.91  63752       44075           44075
## S_262         0.99  42784      0.90  42784       33102           33102
## S_263         0.99  53184      0.90  53184       39177           39177
## S_264         0.99  51337      0.92  51337       36374           36371
## S_265         0.99  47000      0.91  47000       33074           33074
## S_267         0.98  33276      0.88  33276       25918           25918
## S_268         0.98  58828      0.90  58828       44286           44286
## S_269         0.99  55474      0.92  55474       41040           41040
## S_270         0.99  47109      0.93  47109       33787           33787
## S_271         0.99  62749      0.94  62749       49311           49311
## S_272         0.99  63342      0.93  63342       47625           47625
## S_273         0.99  48976      0.94  48976       38035           38035
## S_274         0.99  80614      0.93  80614       59937           59937
## S_275         0.99  42485      0.92  42485       30890           30890
## S_276         0.99  47698      0.92  47698       36986           36986
## S_277         0.98  52627      0.92  52627       40691           40691
## S_278         0.98  41374      0.91  41374       31897           31897
## S_279         0.99  64030      0.93  64030       47728           47728
## S_280         0.99  47420      0.93  47420       37306           37306
## S_281         0.99  43385      0.93  43385       32803           32803
## S_298         0.99  41344      0.93  41344       32872           32872
## S_299         0.99  32859      0.91  32859       25123           25123
## S_300         0.99  44429      0.93  44429       34215           34215
## S_301         0.99  50825      0.93  50825       39342           39342
## S_302         0.99  36206      0.94  36206       27584           27584
## S_303         0.99  51798      0.94  51798       36753           36753
## S_304         0.99  50705      0.93  50705       34090           34090
## S_305         0.99  43798      0.91  43798       32369           32369
## S_306         0.99  46319      0.92  46319       33486           33486
## S_308         0.98  51026      0.85  51026       37026           37026
## S_309         0.99  43102      0.89  43102       30284           30284
## S_310         0.99  33431      0.86  33431       24770           24770
## S_311         0.99  58799      0.88  58799       41706           41706
## S_312         0.99  47333      0.88  47333       34794           34794
## S_313         0.99  48093      0.94  48093       38160           38160
## S_315         0.99  50234      0.91  50234       37273           37273
## S_316         0.99  48036      0.93  48036       35747           35747
## S_317         0.99  50207      0.92  50207       36364           36364
## S_318         0.99  43142      0.93  43142       32689           32689
## S_319         0.99  46711      0.93  46711       35016           35016
## S_320         0.99  45769      0.92  45769       33891           33891
## S_321         0.99  48984      0.92  48984       35926           35926
## S_322         0.99  47055      0.94  47055       33868           33868
## S_323         0.99  50335      0.92  50335       38224           38224
## S_324         0.99  49621      0.95  49621       41368           41368
## S_325         0.99  39705      0.93  39705       29213           29213
## S_326         0.99  49886      0.92  49886       37943           37943
## S_327         0.99  41415      0.94  41415       33353           33353
## S_341         0.99  47667      0.94  47667       37388           37388
## S_342         0.99  36114      0.95  36114       30253           30253
## S_343         0.99  28431      0.92  28431       22578           22578
## S_344         0.99  52379      0.93  52379       41121           41121
## S_345         0.99  23818      0.90  23818       18120           18120
## S_346         0.99  38771      0.92  38771       30174           30174
## S_347         0.99  45983      0.94  45983       34059           34059
## S_348         0.99  46988      0.96  46988       36536           36536
## S_349         0.99  75663      0.94  75663       54422           54422
## S_350         0.99  48047      0.93  48047       33745           33745
## S_351         0.99  67645      0.93  67645       47662           47662
## S_352         0.99  63239      0.92  63239       45958           45958
## S_353         0.99  56477      0.93  56477       45184           45184
## S_354         0.99  48051      0.94  48051       34946           34946
## S_355         0.99  46169      0.92  46169       32861           32861
## S_356         0.99  41514      0.90  41514       29841           29841
## S_438         0.99  60561      0.93  60561       52584           52584
## S_439         0.99  60150      0.93  60150       50751           50751
## S_440         0.99  38560      0.93  38560       33283           33283
## S_441         0.99  66477      0.95  66477       58699           58699
## S_442         0.99  60377      0.96  60377       58881           58881
## S_443         0.99  61479      0.94  61479       50758           50758
## S_444         0.99  63124      0.96  63124       60863           60863
## S_445         0.99  42723      0.94  42723       36375           36375
## S_446         0.99  61694      0.94  61694       50918           50918
## S_447         0.99  63155      0.96  63155       58042           58042
## S_448         0.99  54043      0.95  54043       46975           46975
## S_449         0.99  65940      0.96  65940       57531           57531
## S_450         0.99  62165      0.95  62165       53894           53894
## S_465         0.99  67008      0.97  67008       63618           63618
## S_466         0.99  68679      0.95  68679       60000           60000
## S_467         0.99  59847      0.94  59847       49010           49010
## S_468         0.99  63359      0.93  63359       53096           53096
## S_469         0.99  19263      0.89  19263       15459           15459
## S_470         0.99  65289      0.95  65289       57063           57063
## S_471         0.99  58034      0.89  58034       43164           43164
## S_472         0.99  57373      0.91  57373       43837           43837
## S_473         0.99  72195      0.95  72195       56134           56134
## S_474         0.99  52673      0.88  52673       37386           37384
## S_475         0.99  60887      0.95  60887       48463           48463
## S_476         0.99  55620      0.91  55620       43056           43056
## S_477         0.99  35722      0.91  35722       27699           27699
## S_478         0.99  57328      0.86  57328       43663           43663
##       tabled_pc chimera_out_pc length_filtered_pc    day treatment mouse_label
## S_219         1           0.71                  1 Day_-1    TreatA         m_8
## S_220         1           0.69                  1 Day_-1    TreatA        m_19
## S_221         1           0.73                  1 Day_-1    TreatA        m_34
## S_222         1           0.77                  1 Day_-1    TreatA        m_26
## S_223         1           0.73                  1 Day_-1    TreatA        m_39
## S_224         1           0.73                  1 Day_-1    TreatA        m_46
## S_225         1           0.72                  1 Day_-1    TreatA        m_48
## S_226         1           0.73                  1 Day_-1    TreatA        m_56
## S_227         1           0.72                  1 Day_-1    TreatB        m_29
## S_228         1           0.75                  1 Day_-1    TreatB        m_13
## S_229         1           0.83                  1 Day_-1    TreatB        m_18
## S_230         1           0.73                  1 Day_-1    TreatB        m_27
## S_231         1           0.79                  1 Day_-1    TreatB        m_32
## S_232         1           0.70                  1 Day_-1    TreatB        m_49
## S_233         1           0.76                  1 Day_-1    TreatB        m_54
## S_234         1           0.72                  1 Day_-1    TreatB        m_35
## S_251         1           0.85                  1 Day_-1    TreatC         m_7
## S_252         1           0.74                  1 Day_-1    TreatC        m_10
## S_253         1           0.70                  1 Day_-1    TreatC        m_11
## S_254         1           0.72                  1 Day_-1    TreatC        m_24
## S_255         1           0.74                  1 Day_-1    TreatC        m_28
## S_256         1           0.75                  1 Day_-1    TreatC        m_36
## S_257         1           0.74                  1 Day_-1    TreatC        m_51
## S_259         1           0.76                  1 Day_-1   Control         m_2
## S_260         1           0.70                  1 Day_-1   Control        m_21
## S_261         1           0.69                  1 Day_-1   Control        m_22
## S_262         1           0.77                  1 Day_-1   Control        m_25
## S_263         1           0.74                  1 Day_-1   Control        m_37
## S_264         1           0.71                  1 Day_-1   Control        m_44
## S_265         1           0.70                  1 Day_-1   Control        m_52
## S_267         1           0.78                  1 Day_12    TreatA         m_8
## S_268         1           0.75                  1 Day_12    TreatA        m_19
## S_269         1           0.74                  1 Day_12    TreatA        m_34
## S_270         1           0.72                  1 Day_12    TreatA        m_26
## S_271         1           0.79                  1 Day_12    TreatA        m_39
## S_272         1           0.75                  1 Day_12    TreatA        m_46
## S_273         1           0.78                  1 Day_12    TreatA        m_48
## S_274         1           0.74                  1 Day_12    TreatB        m_29
## S_275         1           0.73                  1 Day_12    TreatB        m_13
## S_276         1           0.78                  1 Day_12    TreatB        m_18
## S_277         1           0.77                  1 Day_12    TreatB        m_27
## S_278         1           0.77                  1 Day_12    TreatB        m_32
## S_279         1           0.75                  1 Day_12    TreatB        m_49
## S_280         1           0.79                  1 Day_12    TreatB        m_54
## S_281         1           0.76                  1 Day_12    TreatB        m_35
## S_298         1           0.80                  1 Day_12    TreatC         m_7
## S_299         1           0.76                  1 Day_12    TreatC        m_10
## S_300         1           0.77                  1 Day_12    TreatC        m_11
## S_301         1           0.77                  1 Day_12    TreatC        m_24
## S_302         1           0.76                  1 Day_12    TreatC        m_28
## S_303         1           0.71                  1 Day_12    TreatC        m_36
## S_304         1           0.67                  1 Day_12    TreatC        m_51
## S_305         1           0.74                  1 Day_12   Control         m_2
## S_306         1           0.72                  1 Day_12   Control        m_21
## S_308         1           0.73                  1 Day_12   Control        m_25
## S_309         1           0.70                  1 Day_12   Control        m_37
## S_310         1           0.74                  1 Day_12   Control        m_44
## S_311         1           0.71                  1 Day_12   Control        m_52
## S_312         1           0.74                  1 Day_12   Control        m_53
## S_313         1           0.79                  1 Day_14    TreatA        m_19
## S_315         1           0.74                  1 Day_14    TreatA        m_34
## S_316         1           0.74                  1 Day_14    TreatA        m_26
## S_317         1           0.72                  1 Day_14    TreatA        m_39
## S_318         1           0.76                  1 Day_14    TreatA        m_46
## S_319         1           0.75                  1 Day_14    TreatA        m_48
## S_320         1           0.74                  1 Day_14    TreatB        m_29
## S_321         1           0.73                  1 Day_14    TreatB        m_13
## S_322         1           0.72                  1 Day_14    TreatB        m_18
## S_323         1           0.76                  1 Day_14    TreatB        m_27
## S_324         1           0.83                  1 Day_14    TreatB        m_32
## S_325         1           0.74                  1 Day_14    TreatB        m_49
## S_326         1           0.76                  1 Day_14    TreatB        m_54
## S_327         1           0.81                  1 Day_14    TreatB        m_35
## S_341         1           0.78                  1 Day_14    TreatC         m_7
## S_342         1           0.84                  1 Day_14    TreatC        m_10
## S_343         1           0.79                  1 Day_14    TreatC        m_11
## S_344         1           0.79                  1 Day_14    TreatC        m_24
## S_345         1           0.76                  1 Day_14    TreatC        m_28
## S_346         1           0.78                  1 Day_14    TreatC        m_36
## S_347         1           0.74                  1 Day_14    TreatC        m_51
## S_348         1           0.78                  1 Day_14    TreatC        m_55
## S_349         1           0.72                  1 Day_14   Control         m_2
## S_350         1           0.70                  1 Day_14   Control        m_21
## S_351         1           0.70                  1 Day_14   Control        m_22
## S_352         1           0.73                  1 Day_14   Control        m_25
## S_353         1           0.80                  1 Day_14   Control        m_37
## S_354         1           0.73                  1 Day_14   Control        m_44
## S_355         1           0.71                  1 Day_14   Control        m_52
## S_356         1           0.72                  1 Day_14   Control        m_53
## S_438         1           0.87                  1  Day_6    TreatA        m_19
## S_439         1           0.84                  1  Day_6    TreatA        m_26
## S_440         1           0.86                  1  Day_6    TreatA        m_46
## S_441         1           0.88                  1  Day_6    TreatA        m_48
## S_442         1           0.98                  1  Day_6    TreatA        m_56
## S_443         1           0.83                  1  Day_6    TreatB        m_29
## S_444         1           0.96                  1  Day_6    TreatB        m_13
## S_445         1           0.85                  1  Day_6    TreatB        m_18
## S_446         1           0.83                  1  Day_6    TreatB        m_27
## S_447         1           0.92                  1  Day_6    TreatB        m_32
## S_448         1           0.87                  1  Day_6    TreatB        m_49
## S_449         1           0.87                  1  Day_6    TreatB        m_54
## S_450         1           0.87                  1  Day_6    TreatB        m_35
## S_465         1           0.95                  1  Day_6    TreatC         m_7
## S_466         1           0.87                  1  Day_6    TreatC        m_10
## S_467         1           0.82                  1  Day_6    TreatC        m_11
## S_468         1           0.84                  1  Day_6    TreatC        m_24
## S_469         1           0.80                  1  Day_6    TreatC        m_28
## S_470         1           0.87                  1  Day_6    TreatC        m_36
## S_471         1           0.74                  1  Day_6   Control         m_2
## S_472         1           0.76                  1  Day_6   Control        m_21
## S_473         1           0.78                  1  Day_6   Control        m_22
## S_474         1           0.71                  1  Day_6   Control        m_25
## S_475         1           0.80                  1  Day_6   Control        m_37
## S_476         1           0.77                  1  Day_6   Control        m_44
## S_477         1           0.78                  1  Day_6   Control        m_52
## S_478         1           0.76                  1  Day_6   Control        m_53
##          BW   BW_percent BW_delta DAI treatment_grouped treatment_invivo
## S_219 20.55  0.000000000     0.00   0            TreatA             none
## S_220 21.32  0.000000000     0.00   0            TreatA             none
## S_221 20.11  0.000000000     0.00   0            TreatA             none
## S_222 21.90  0.000000000     0.00   0            TreatA             none
## S_223 20.94  0.000000000     0.00   0            TreatA             none
## S_224 21.12  0.000000000     0.00   0            TreatA             none
## S_225 20.35  0.000000000     0.00   0            TreatA             none
## S_226 20.51  0.000000000     0.00   0            TreatA             none
## S_227 20.95  0.000000000     0.00   0            TreatB             none
## S_228 20.35  0.000000000     0.00   0            TreatB             none
## S_229 21.55  0.000000000     0.00   0            TreatB             none
## S_230 23.18  0.000000000     0.00   0            TreatB             none
## S_231 20.78  0.000000000     0.00   0            TreatB             none
## S_232 20.90  0.000000000     0.00   0            TreatB             none
## S_233 19.70  0.000000000     0.00   0            TreatB             none
## S_234 20.53  0.000000000     0.00   0            TreatB             none
## S_251 20.98  0.000000000     0.00   0            TreatC             none
## S_252 20.48  0.000000000     0.00   0            TreatC             none
## S_253 20.47  0.000000000     0.00   0            TreatC             none
## S_254 22.91  0.000000000     0.00   0            TreatC             none
## S_255 21.30  0.000000000     0.00   0            TreatC             none
## S_256 22.47  0.000000000     0.00   0            TreatC             none
## S_257 21.92  0.000000000     0.00   0            TreatC             none
## S_259 21.09  0.000000000     0.00   0           Control             none
## S_260 21.31  0.000000000     0.00   0           Control             none
## S_261 20.45  0.000000000     0.00   0           Control             none
## S_262 22.00  0.000000000     0.00   0           Control             none
## S_263 20.38  0.000000000     0.00   0           Control             none
## S_264 21.60  0.000000000     0.00   0           Control             none
## S_265 21.75  0.000000000     0.00   0           Control             none
## S_267 18.24  0.112408759     0.20   1            TreatA           TreatA
## S_268 18.64  0.125703565     0.64   1            TreatA           TreatA
## S_269 19.75  0.017901542     0.90   0            TreatA           TreatA
## S_270 22.45 -0.025114155    -0.15   1            TreatA           TreatA
## S_271 20.86  0.003820439     1.81   0            TreatA           TreatA
## S_272 20.12  0.047348485     0.70   0            TreatA           TreatA
## S_273 19.75  0.029484029     0.71   0            TreatA           TreatA
## S_274 17.62  0.158949881     0.26   3            TreatB           TreatB
## S_275 17.60  0.135135135     0.05   3            TreatB           TreatB
## S_276 17.57  0.184686775     0.39   3            TreatB           TreatB
## S_277 20.47  0.116911130    -0.46   1            TreatB           TreatB
## S_278 17.32  0.166506256     0.24   3            TreatB           TreatB
## S_279 17.95  0.141148325     0.18   3            TreatB           TreatB
## S_280 17.25  0.124365482     0.20   3            TreatB           TreatB
## S_281 18.45  0.101315149    -0.10   3            TreatB           TreatB
## S_298 18.25  0.130123928     0.30   3            TreatC           TreatC
## S_299 17.35  0.152832031     0.25   3            TreatC           TreatC
## S_300 16.85  0.176844162    -0.01   4            TreatC           TreatC
## S_301 20.01  0.126582278    -0.35   2            TreatC           TreatC
## S_302 18.12  0.149295775    -0.05   3            TreatC           TreatC
## S_303 17.07  0.240320427    -0.19   4            TreatC           TreatC
## S_304 16.07  0.266879562     0.20   4            TreatC           TreatC
## S_305 20.70 -0.018492176    -0.02   0           Control             none
## S_306 21.54  0.010793055     0.74   0           Control             none
## S_308 21.86 -0.006363636    -1.88   0           Control             none
## S_309 21.16  0.038272816    -0.22   0           Control             none
## S_310 22.97  0.063425926     0.12   0           Control             none
## S_311 22.07  0.014712644     0.67   0           Control             none
## S_312 21.47 -0.029823769    -0.87   0           Control             none
## S_313 21.13  0.008911820     1.32   0            TreatA           TreatA
## S_315 21.03 -0.045748384     0.84   0            TreatA           TreatA
## S_316 21.34  0.025570776    -0.85   0            TreatA           TreatA
## S_317 22.88 -0.092645654     2.42   0            TreatA           TreatA
## S_318 21.81 -0.032670455     0.59   0            TreatA           TreatA
## S_319 21.31 -0.047174447     0.49   0            TreatA           TreatA
## S_320 18.02  0.139856802     0.17   2            TreatB           TreatB
## S_321 18.23  0.104176904     0.38   2            TreatB           TreatB
## S_322 18.01  0.164269142     0.45   3            TreatB           TreatB
## S_323 20.95  0.096203624     0.94   0            TreatB           TreatB
## S_324 17.85  0.141000962     0.50   2            TreatB           TreatB
## S_325 18.35  0.122009569     0.23   2            TreatB           TreatB
## S_326 18.52  0.059898477     0.52   0            TreatB           TreatB
## S_327 18.95  0.076960546     0.41   0            TreatB           TreatB
## S_341 18.21  0.132030505     0.07   2            TreatC           TreatC
## S_342 17.58  0.141601563    -0.07   3            TreatC           TreatC
## S_343 17.65  0.137762579     0.09   2            TreatC           TreatC
## S_344 20.23  0.116979485    -0.40   3            TreatC           TreatC
## S_345 18.96  0.109859155     0.42   3            TreatC           TreatC
## S_346 18.02  0.198041834     0.80   3            TreatC           TreatC
## S_347 16.96  0.226277372     0.72   4            TreatC           TreatC
## S_348 18.02  0.121404193     0.06   2            TreatC           TreatC
## S_349 21.28  0.009009009     0.58   0           Control             none
## S_350 21.72  0.019239794     0.18   0           Control             none
## S_351 21.70  0.061124694     0.54   0           Control             none
## S_352 21.79 -0.009545455    -0.07   0           Control             none
## S_353 21.72  0.065750736     0.56   0           Control             none
## S_354 23.77  0.100462963     0.80   0           Control             none
## S_355 21.75  0.000000000    -0.26   0           Control             none
## S_356 22.26  0.005874379     0.79   0           Control             none
## S_438 20.48  0.039399625    -0.93   1            TreatA           TreatA
## S_439 20.35  0.070776256    -1.61   2            TreatA           TreatA
## S_440 20.77  0.016571970    -0.49   1            TreatA           TreatA
## S_441 19.90  0.022113022    -0.43   1            TreatA           TreatA
## S_442 19.50  0.049244271    -1.12   1            TreatA           TreatA
## S_443 21.32 -0.017661098     0.39   2            TreatB           TreatB
## S_444 20.59 -0.011793612    -0.79   1            TreatB           TreatB
## S_445 21.06  0.022737819    -1.04   1            TreatB           TreatB
## S_446 22.95  0.009922347    -0.21   0            TreatB           TreatB
## S_447 20.91 -0.006256015    -0.15   1            TreatB           TreatB
## S_448 20.41  0.023444976    -1.00   2            TreatB           TreatB
## S_449 19.93 -0.011675127    -0.86   1            TreatB           TreatB
## S_450 20.81 -0.013638578    -0.60   1            TreatB           TreatB
## S_465 20.98  0.000000000    -0.93   1            TreatC           TreatC
## S_466 20.42  0.002929687    -0.69   1            TreatC           TreatC
## S_467 19.93  0.026380068    -0.87   1            TreatC           TreatC
## S_468 23.81 -0.039284155     0.46   0            TreatC           TreatC
## S_469 20.73  0.026760563    -0.70   1            TreatC           TreatC
## S_470 21.79  0.030262572    -0.33   2            TreatC           TreatC
## S_471 20.65 -0.020862968    -0.34   0           Control             none
## S_472 21.20 -0.005161896     0.25   0           Control             none
## S_473 20.90  0.022004890    -1.05   0           Control             none
## S_474 21.68 -0.014545455     0.38   0           Control             none
## S_475 21.68  0.063788027     0.43   0           Control             none
## S_476 22.45  0.039351852     0.52   0           Control             none
## S_477 20.76 -0.045517241    -0.19   0           Control             none
## S_478 21.74 -0.017623136     0.18   0           Control             none
##       day_num
## S_219      -1
## S_220      -1
## S_221      -1
## S_222      -1
## S_223      -1
## S_224      -1
## S_225      -1
## S_226      -1
## S_227      -1
## S_228      -1
## S_229      -1
## S_230      -1
## S_231      -1
## S_232      -1
## S_233      -1
## S_234      -1
## S_251      -1
## S_252      -1
## S_253      -1
## S_254      -1
## S_255      -1
## S_256      -1
## S_257      -1
## S_259      -1
## S_260      -1
## S_261      -1
## S_262      -1
## S_263      -1
## S_264      -1
## S_265      -1
## S_267      12
## S_268      12
## S_269      12
## S_270      12
## S_271      12
## S_272      12
## S_273      12
## S_274      12
## S_275      12
## S_276      12
## S_277      12
## S_278      12
## S_279      12
## S_280      12
## S_281      12
## S_298      12
## S_299      12
## S_300      12
## S_301      12
## S_302      12
## S_303      12
## S_304      12
## S_305      12
## S_306      12
## S_308      12
## S_309      12
## S_310      12
## S_311      12
## S_312      12
## S_313      14
## S_315      14
## S_316      14
## S_317      14
## S_318      14
## S_319      14
## S_320      14
## S_321      14
## S_322      14
## S_323      14
## S_324      14
## S_325      14
## S_326      14
## S_327      14
## S_341      14
## S_342      14
## S_343      14
## S_344      14
## S_345      14
## S_346      14
## S_347      14
## S_348      14
## S_349      14
## S_350      14
## S_351      14
## S_352      14
## S_353      14
## S_354      14
## S_355      14
## S_356      14
## S_438       6
## S_439       6
## S_440       6
## S_441       6
## S_442       6
## S_443       6
## S_444       6
## S_445       6
## S_446       6
## S_447       6
## S_448       6
## S_449       6
## S_450       6
## S_465       6
## S_466       6
## S_467       6
## S_468       6
## S_469       6
## S_470       6
## S_471       6
## S_472       6
## S_473       6
## S_474       6
## S_475       6
## S_476       6
## S_477       6
## S_478       6
```

```r
ps_invivo %>% 
    filter_taxa(function(x) sum(x > 0) > 0, TRUE) %>% 
  saveRDS("~/Documents/GitHub/DivComAnalyses/data-raw/ps_invivo.RDS")
```

```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] phyloseq_1.34.0 forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7    
##  [5] purrr_0.3.4     readr_2.0.1     tidyr_1.1.4     tibble_3.1.6   
##  [9] ggplot2_3.3.5   tidyverse_1.3.1
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-152        fs_1.5.2            lubridate_1.8.0    
##  [4] httr_1.4.2          tools_4.0.2         backports_1.4.1    
##  [7] bslib_0.2.5.1       vegan_2.5-7         utf8_1.2.2         
## [10] R6_2.5.1            mgcv_1.8-36         DBI_1.1.1          
## [13] BiocGenerics_0.34.0 colorspace_2.0-2    permute_0.9-5      
## [16] ade4_1.7-17         withr_2.4.3         tidyselect_1.1.1   
## [19] compiler_4.0.2      cli_3.1.0           rvest_1.0.1        
## [22] Biobase_2.50.0      xml2_1.3.2          sass_0.4.1         
## [25] scales_1.1.1        digest_0.6.29       rmarkdown_2.13     
## [28] XVector_0.28.0      pkgconfig_2.0.3     htmltools_0.5.2    
## [31] dbplyr_2.1.1        fastmap_1.1.0       rlang_1.0.2        
## [34] readxl_1.3.1        rstudioapi_0.13     jquerylib_0.1.4    
## [37] generics_0.1.1      jsonlite_1.8.0      magrittr_2.0.3     
## [40] biomformat_1.16.0   Matrix_1.3-4        Rcpp_1.0.7         
## [43] munsell_0.5.0       S4Vectors_0.26.1    Rhdf5lib_1.10.1    
## [46] fansi_0.5.0         ape_5.5             lifecycle_1.0.1    
## [49] stringi_1.7.6       yaml_2.3.5          MASS_7.3-54        
## [52] zlibbioc_1.34.0     rhdf5_2.32.4        plyr_1.8.6         
## [55] grid_4.0.2          parallel_4.0.2      crayon_1.4.2       
## [58] lattice_0.20-44     splines_4.0.2       Biostrings_2.56.0  
## [61] haven_2.4.3         multtest_2.44.0     hms_1.1.0          
## [64] knitr_1.38          pillar_1.6.4        igraph_1.2.6       
## [67] reshape2_1.4.4      codetools_0.2-18    stats4_4.0.2       
## [70] reprex_2.0.1        glue_1.6.2          evaluate_0.14      
## [73] data.table_1.14.2   modelr_0.1.8        vctrs_0.3.8        
## [76] tzdb_0.1.2          foreach_1.5.1       cellranger_1.1.0   
## [79] gtable_0.3.0        assertthat_0.2.1    xfun_0.30          
## [82] broom_0.7.11        survival_3.2-13     iterators_1.0.13   
## [85] IRanges_2.22.2      cluster_2.1.2       ellipsis_0.3.2
```

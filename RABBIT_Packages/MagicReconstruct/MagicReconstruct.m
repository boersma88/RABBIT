(* Mathematica Package *)

(* Created by the Wolfram Workbench Mar 18, 2014 *)

BeginPackage["MagicReconstruct`",{"ContinuousTimeHmm`","MagicReconstruct`DataValidation`","MagicModel`","MagicAlgorithm`"}]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["MagicReconstruct`*"];
ClearAll @@ Names["MagicReconstruct`*"];

HMMMethod::usage = "HMMMethod is an option to specify the alogrithm used for Hidden Markov Model, it has to be one of \"origPathSampling\", \"origPosteriorDecoding\", or \"origViterbiDecoding\". "

SampleSize::usage = "SampleSize is an option to specify the sample size when by default Method->origPathSampling. By default, SampleSize->1000."

PrintTimeElapsed::usage = "PrintTimeElapsed is an option to secify whether to print time elapsed."

magicReconstruct::usage = "magicReconstruct[magicSNP, model, epsF, eps, matingScheme, outputfile] reconstructs the ancestral origins in multi-parental populations. magicSNP is the input marker data, which are analyzed under the model using a HMM algorithm HMMMethod. The parameters epsF and eps refer to the allelic typing errors for founders and sampled individuals, respectively. The parameter matingScheme is a list of mating schemes of the mapping population until the last generation. The results are saved in outputfile. \n magicReconstruct[magicSNPfile, model, epsF, eps, matingScheme, outfile] where the input marker data are privided as the filename magicSNPfile."

saveAsSummaryMR::usage = "saveAsSummaryMR[resultFile, summaryFile] produces a summary of the output results from magicReconstruct, and save in summaryFile."

getSummaryMR::usage = "getSummaryMR[summaryFile] imports the summaryFile produced by saveAsSummaryMR."

origGenotypes::usage = "origGenotypes[nFounder] gives {{genotypes, diplotypes}, {geno2diplo, diplo2geno}}.  The genotypes and diplotypes are the definitions where the inbred founders are labeled by natural numbers starting from 1,geno2diplo gives the mapping from genotypes to dipotypes, and diplo2geno gives the mapping from diplotypes to genotypes."

bestAncestry::usage = "bestAncestry[prob] gives the ancestral state calls by maximum posterior diplotype/genotype/haplotype probabilities.  "

toGenoprob::usage = "toGenoprob[diploprob] transfers conditional probabilities from L^2 diplotypes into L(L+1)/2 genotypes, where L is the number of inbred founder."

toHaploprob::usage = "toHaploprob[genoprob] transfers conditional probabilities from L(L+1)/2 genotypes into L haplotypes, where L is the number of inbred founder."

calOrigLogl::usage = "calOrigLogl[magicSNP, model, epsF, eps, matingScheme] calculates the marginal likelihood. "

calOrigGeneration::usage = "calOrigGeneration[magicSNP, model, epsF, eps, Fmin, matingSchemeMax] calculates the posterior probabilities of the sample generation. "

(*transformMagicSNP::usage ="transformMagicSNP"*)

Begin["`Private`"]
(* Implementation of the package *)

(*diplotypes refer to a single locus, phased genotypes*)
(*geno2diplo[[i]] refer to the index of the ith genotype, in the diplotypes*)
origGenotypes[nFounder_] :=
    Module[ {genotypes,diplotypes,geno2diplo,diplo2geno, i,j},
        diplotypes = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1];
        (*lower triangular matrix*)
        genotypes = Flatten[Table[{i, j}, {i, nFounder}, {j, i}], 1];
        diplo2geno = Flatten[Position[genotypes, #, {1}, 1] & /@ (Sort[#, Greater] & /@diplotypes)];
        geno2diplo = {#[[All, 1]], #[[1, 2]]} & /@ 
            SplitBy[SortBy[Transpose[{Range[Length[diplotypes]], diplo2geno}], Last],Last];
        geno2diplo = geno2diplo[[All, 1]];
        {{genotypes,diplotypes},{geno2diplo,diplo2geno}}
    ]

(*diploprob[[i,j]]: a list of diplotype probabilities at marker j of linkage group i, for a given indvidual*)
(*also works if there are no level of linkage groups, also works if there are additonal level of individuals*)    
toGenoprob[diploprob_] :=
    Module[ {depth,len,nFounder, geno2diplo, genoprob, f},
        depth = Depth[diploprob];
        len = Union[Flatten[Map[Length, diploprob, {depth - 2}]]];
        nFounder = Sqrt[First[len]];
        If[ ! (Length[len] == 1 && IntegerQ[nFounder]),
            Print["toGenoprob: the argument is not an array of diplotype probabilities!"];
            Abort[]
        ];
        geno2diplo = origGenotypes[nFounder][[2,1]];
        f = Function[{x}, Total[x[[#]] & /@ geno2diplo, {2}]];
        genoprob = Map[f, diploprob, {depth - 2}];
        genoprob
    ]

toHaploprob[genoprob_] :=
    Module[ {depth, len, nFounder, genotypes, mtx,haploprob,x,i,j},
        depth = Depth[genoprob];
        len = Union[Flatten[Map[Length, genoprob, {depth - 2}]]];
        nFounder = Solve[x (x + 1)/2 == First[len] && x > 0, x][[1, 1, 2]];
        If[ ! (Length[len] == 1 && IntegerQ[nFounder]),
            Print["toGenoprob: the argument is not an array of unphased genotype probabilities!"];
            Abort[]
        ];
        genotypes = origGenotypes[nFounder][[1, 1]];
        mtx = ConstantArray[0, {Length[genotypes], nFounder}];
        Do[mtx[[i, genotypes[[i, j]]]] += 1, {i, Length[genotypes]}, {j,Dimensions[genotypes][[2]]}];
        mtx = Transpose[mtx]/Dimensions[genotypes][[2]];
        haploprob = Map[mtx.# &, genoprob, {depth - 2}];
        haploprob
    ]

(*prob[[i,j]]: a list of (haplotype/genotype/diplotype) probabilities at marker j of linkage group i for a given individual*)
(*return res[[i,j]] ={indices,max}, the positions of the max probability, for prob[[i,j]]*)
(*also works if there are no level of linkage groups, also works if there are additonal level of individuals*)  
bestAncestry[prob_] :=
    Map[{N[Max[#]], Flatten[Position[#, Max[#], {1}, Heads -> False]]} &, Round[prob, 10^(-5)], {Depth[prob] - 2}]
 

Options[magicReconstruct] = {
    HMMMethod -> "origPathSampling",
    SampleSize -> 1000,
    PrintTimeElapsed ->False
}

magicReconstruct::wrongMethod :=
    "The option HMMMethod has to take \"origPathSampling\", \"origPosteriorDecoding\", or \"origViterbiDecoding\"."

magicReconstruct::wrongModel :=
    "The 2nd model input parameter has to take \"jointModel\", \"indepModel\", or \"depModel\"."
           
magicReconstruct[magicSNPfile_String, model_String, epsF_?NumericQ, eps_?NumericQ,
    matingScheme_?(VectorQ[#, StringQ]&), outfile_String, opts : OptionsPattern[]] :=
    Module[ {magicSNP},
        If[ !MemberQ[{"jointModel","indepModel","depModel"},model],
            Message[magicReconstruct::wrongModel];
            Return[$Failed]
        ];
        magicSNP = Import[magicSNPfile,"CSV"];
        SNPValidation[magicSNP];
        magicReconstruct[magicSNP,model, epsF, eps, matingScheme, outfile, opts]
    ]  

magicReconstruct[magicSNP_List,model_String, epsF_?NumericQ, eps_?NumericQ,
    matingScheme_?(VectorQ[#, StringQ]&), outputfile_String, opts : OptionsPattern[]] :=
    Module[ {algorithmname,samplesize, isprint,starttime,deltd, nFounder,posA,posX,genders,founderid,sampleid,snpMap,outstream,
            startProb, tranProb, startProbAutosome,tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale,
            founderHaplo, obsGeno, dataProb, res,ind,diplotypes,haplotodiplo,prob},
        {algorithmname, samplesize,isprint} = OptionValue@{HMMMethod, SampleSize,PrintTimeElapsed};
        If[ !MemberQ[{"jointModel","indepModel","depModel"},model],
            Message[magicReconstruct::wrongModel];
            Return[$Failed]
        ];
        If[ MemberQ[{"origPathSampling", "origPosteriorDecoding", "origViterbiDecoding"},algorithmname],
            origalgorithm = Symbol[algorithmname],
            Message[magicReconstruct::wrongMethod];
            Return[$Failed]
        ];
        SNPValidation[magicSNP];
        If[ isprint,
            starttime = TimeUsed[];
            Print["Start date =",DateString[]];
        ];
        If[ isprint,
            Print["Pre-computing data probability...","\toutputfile = ",outputfile]
        ];
        {deltd, founderHaplo, obsGeno,genders,posA,posX,nFounder,founderid,sampleid,snpMap} = transformMagicSNP[magicSNP];
        {startProbAutosome, tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale} = 
            magicPriorProcess[snpMap,model,nFounder,matingScheme];
        startProb = tranProb = Table[0,{Length[posA]+Length[posX]}];
        If[ posA=!={},
            {startProb[[posA]], tranProb[[posA]]} = {startProbAutosome, tranProbAutosome};
        ];
        Put[{model,epsF, eps, matingScheme,algorithmname,genders,posA,posX,founderid,sampleid,snpMap}, outputfile];
        Quiet[Close[outputfile]];
        outstream = OpenAppend[outputfile];
        diplotypes = origGenotypes[nFounder][[1, 2]];
        haplotodiplo = Flatten[Position[Equal @@ # & /@ diplotypes, True]];
        Monitor[        
            Do[
             If[ isprint,
                 Print["Time elapsed = ", Round[TimeUsed[] - starttime,0.1], " Seconds.", "\tStarting analyzing line ", ind, " of ", Length[obsGeno]];
             ];
             Switch[genders[[ind]],
                 "Female",
                 {startProb[[posX]], tranProb[[posX]]} = {startProbFemale, tranProbFemale},
                 "Male",
                 {startProb[[posX]], tranProb[[posX]]} = {startProbMale, tranProbMale}                 
             ];
             dataProb = magicSNPLikelihood[founderHaplo, obsGeno[[ind]], epsF, eps, posA, posX, genders[[ind]]];
             res = If[ algorithmname === "origPathSampling",
                       origalgorithm[startProb, tranProb, dataProb, samplesize],
                       origalgorithm[startProb, tranProb, dataProb]
                   ];
             If[ genders[[ind]] == "Male",
                 Switch[
                  algorithmname,
                  "origPosteriorDecoding",
                  prob = ConstantArray[0, ReplacePart[Dimensions[res[[posX, 2]]], -1 -> Length[diplotypes]]];
                  prob[[All, All, haplotodiplo]] = res[[posX, 2]];
                  res[[posX, 2]] = prob,
                  "origViterbiDecoding",
                  res[[posX, 2, ;; -2, 2]] = haplotodiplo[[#]] & /@ res[[posX, 2, ;; -2, 2]],
                  "origPathSampling",
                  res[[posX, 2, All, ;; -2, 2]] = Map[haplotodiplo[[#]] &, res[[posX, 2, All, ;; -2, 2]], {2}],
                  _,
                  Print["Wrong HMMMethod " <> algorithmname <> "!"];
                  ]
             ];
             PutAppend[res, outstream], {ind, Length[obsGeno]}];
             , ProgressIndicator[ind, {1, Length[obsGeno]}]
         ];
        Close[outstream];
        (*summaryfile=StringDrop[outputfile, -4] <> "_Summary.csv";
        If[ isprint,
            Print["Time elapsed = ", Round[TimeUsed[] - starttime,0.1], " Seconds.", "\tSaving as summary in summaryfile = ",summaryfile]
        ];
        saveAsSummaryMR[outputfile, summaryfile];*)
        If[ isprint,
            Print["Done! Finished date =",DateString[], ". \tTime used in total = ", Round[TimeUsed[] - starttime,0.1], " Seconds."]
        ];
    ]

transformGeno[geno_] :=
    Module[ {strgeno = geno,res},
        strgeno[[3 ;;, 2 ;;]] = Map[ToString, strgeno[[3 ;;, 2 ;;]], {2}];
        res = SplitBy[Transpose[strgeno[[2 ;;, 2 ;;]]], First];
        res = Transpose[Transpose[#] & /@ res[[All, All, 2 ;;]]];
        res
    ]
    
transformMagicSNP[magicSNP_List] :=
    Module[ {snpMap, founderData, obsData,founderid,sampleid,deltd,founderHaplo,obsGeno,nFounder,chrs,posX,posA,genders},
        nFounder = magicSNP[[1,2]];
        snpMap = Transpose[magicSNP[[2 ;; 4]]];
        founderData = Join[magicSNP[[2;;3]],magicSNP[[5 ;; 4 + nFounder]]];
        obsData = Join[magicSNP[[2;;3]],magicSNP[[5 + nFounder ;;]]];
        founderid = ToString[#]&/@founderData[[3 ;;, 1]];
        sampleid = ToString[#]&/@obsData[[3 ;;, 1]];
        (*change the genetic distance from centiMorgan into Morgan*)
        deltd = Differences[#] & /@ (SplitBy[snpMap[[2 ;;, 2 ;;]], First][[All, All, 2]]);
        (*change the genetic distance from centiMorgan into Morgan*)
        deltd = 0.01 deltd;
        founderHaplo = transformGeno[founderData];
        obsGeno = transformGeno[obsData];
        (*extract posA and posX*)
        chrs = Split[founderData[[2, 2 ;;]]][[All, 1]];
        posX = Flatten[Position[chrs, "X"]];
        posA = Complement[Range[Length[chrs]], posX];
        (*extract gender from genotype on sex chromsomes*)
        If[ posX === {},
            genders = Table["Monoecious", {Length[obsData]}],
            genders = Union[Flatten[#]] & /@ obsGeno[[All, posX]];
            genders = (Complement[#, {"1", "2", "N"}] === {} & /@ genders);
            genders = genders /. {True -> "Male", False -> "Female"};
        ];
        (*return*)
        {deltd, founderHaplo, obsGeno,genders,posA,posX,nFounder,founderid,sampleid,Transpose[snpMap]}
    ]
        
calOrigLogl[magicSNP_List, model_String, epsF_?NumericQ, 
  eps_?NumericQ, matingScheme_?(VectorQ[#, StringQ] &)] :=
    Module[ {founderHaplo, obsGeno, ind,deltd,nFounder,founderid,sampleid,snpMap,genders,posA,posX,dataProb, logl,
        startProb,tranProb, startProbAutosome, tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale},
        {deltd, founderHaplo, obsGeno,genders,posA,posX,nFounder,founderid,sampleid,snpMap} = transformMagicSNP[magicSNP];
        {startProbAutosome, tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale} = 
            magicPriorProcess[snpMap,model,nFounder,matingScheme];
        startProb = tranProb = Table[0,{Length[posA]+Length[posX]}];
        If[ posA=!={},
            {startProb[[posA]], tranProb[[posA]]} = {startProbAutosome, tranProbAutosome};
        ];
        logl = Table[
          dataProb = magicSNPLikelihood[founderHaplo, obsGeno[[ind]], epsF, eps, posA, posX, genders[[ind]]];
          Switch[genders[[ind]],
                 "Female",
                 {startProb[[posX]], tranProb[[posX]]} = {startProbFemale, tranProbFemale},
                 "Male",
                 {startProb[[posX]], tranProb[[posX]]} = {startProbMale, tranProbMale}                 
             ];
          origLogLiklihood[startProb, tranProb, dataProb], {ind,Length[obsGeno]}];
        logl
    ]

calOrigGeneration[magicSNP_List, model_String, epsF_?NumericQ, 
  eps_?NumericQ, Fmin_Integer, matingSchemeMax_?(VectorQ[#, StringQ] &)] :=
    Module[ {Fmax, schemelist,i,lls},
        Fmax = Length[matingSchemeMax];
        schemelist = Table[Take[matingSchemeMax, i], {i, Fmin, Fmax}];
        lls = calOrigGeneration[magicSNP, model, epsF, eps, schemelist];
        lls[[All, 2, 1]] = "Generation";
        lls[[All, 3 ;;, 1]] += Fmin - 1;
        lls
    ]

calOrigGeneration[magicSNP_List, model_String, epsF_?NumericQ, 
  eps_?NumericQ, schemelist_?(VectorQ[#, ListQ] &),isSameGeneration_:False] :=
    Module[ {loglls, lls, nFounder, id, i},
        Monitor[
              loglls = Table[calOrigLogl[magicSNP, model, epsF, eps, schemelist[[i]]], {i, Length[schemelist]}], 
                 ProgressIndicator[i, {1, Length[schemelist]}]];
        nFounder = magicSNP[[1, 2]];
        id = magicSNP[[nFounder + 5 ;;, 1]];
        If[ isSameGeneration,
            lls = Total[loglls, {2}];
            lls = Normalize[Exp[lls - Max[lls]], Total];
            lls = Sort[
              Transpose[{Range[Length[schemelist]], lls}], #1[[2]] > #2[[2]] &];
            lls = Join[{{"SchemeIndex", "PosteriorProb"}}, lls],
            lls = Normalize[Exp[# - Max[#]], Total] & /@ Transpose[loglls];
            lls = Sort[Transpose[{Range[Length[schemelist]], #}], #1[[2]] > #2[[2]] &] & /@ lls;
            lls = Transpose[Join[{id, Table[{"SchemeIndex", "PosteriorProb"}, {Length[id]}]}, Transpose[lls]]]
        ]
    ]
  
saveAsSummaryMR[resultFile_String?FileExistsQ, summaryFile_String] :=
    Module[ {res, model, epsF, eps, matingScheme,genders, posA, posX, nFounder, algorithmname, chrid,logl,logl2,
      founderid, sampleid, snpmap, genotypes, diplotypes, genoID, genotypes2, haploID, haplotypes2, 
      genoprob, haploprob, rowID, genoprob2, haploprob2, diploID, diplotypes2, path, summary, diplotohaplo,key = "MagicReconstruct-Summary"},
        res = ReadList[resultFile];
        {model,epsF, eps, matingScheme,algorithmname,genders,posA,posX,founderid,sampleid,snpmap} = res[[1]];
        nFounder = Length[founderid];
        {genotypes, diplotypes} = origGenotypes[nFounder][[1]];
        diplotohaplo = Replace[diplotypes, {{x_, x_} :> x, {_, _} -> 0}, {1}];
        haploID = "haplotype" <> ToString[#] & /@ Range[nFounder];
        haplotypes2 = Transpose[{haploID, Range[nFounder], founderid}];
        haplotypes2 = Join[{{"Haplotype", "Code", "founder"}}, haplotypes2];
        genoID = "genotype" <> ToString[#] & /@ Range[Length[genotypes]];
        genotypes2 = Transpose[{genoID,
           StringJoin[#[[1]], "---", #[[2]]] & /@ Map[ToString, genotypes, {2}],
           StringJoin[#[[1]], "---", #[[2]]] & /@ Map[founderid[[#]] &, genotypes, {2}]}];
        genotypes2 = Join[{{"Genotype", "Code", "founder"}}, genotypes2];
        diploID = "diplotype" <> ToString[#] & /@ Range[Length[diplotypes]];
        diplotypes2 = Transpose[{diploID,
           StringJoin[#[[1]], "---", #[[2]]] & /@ Map[ToString, diplotypes, {2}],
           StringJoin[#[[1]], "---", #[[2]]] & /@ Map[founderid[[#]] &, diplotypes, {2}]}];
        diplotypes2 = Join[{{"Diplotype", "Code", "founder"}}, diplotypes2];
        chrid = Split[snpmap[[2, 2 ;;]]][[All, 1]];
        logl = res[[2 ;;, All, 1]];
        logl2 = Join[Transpose[{Prepend[sampleid, "Lines/Chr"]}],Prepend[logl, chrid], 2];
        Switch[
         algorithmname, 
         "origPosteriorDecoding",
         genoprob = toGenoprob[res[[2 ;;, All, 2]]];
         haploprob = toHaploprob[genoprob];
         (*to output haplotypes2 and genotypes2,*)
         (*to output genoprob2 and haploprob2*)
         If[ model === "depModel",
             genoprob2 = {None},
             rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, genoID]];
             genoprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ genoprob, 1];
             genoprob2 = Join[snpmap,Join[Transpose[{rowID}], genoprob2, 2]];
         ];
         rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, haploID]];
         haploprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ haploprob, 1];
         haploprob2 = Join[snpmap,Join[Transpose[{rowID}], haploprob2, 2]];         
         (*export*)
         summary = Join[{{key, "Genetic map of biallelic markers"}}, snpmap,
           {{key, "Ln marginal likelihood"}}, logl2,
           {{key, "Genotypes in order"}}, genotypes2,
           {{key, "Conditonal genotype probability"}}, genoprob2,
           {{key, "haplotypes in order"}}, haplotypes2,
           {{key, "Conditonal haplotype probability"}}, haploprob2];
         Export[summaryFile, ExportString[summary, "CSV"], "Table"],
         "origViterbiDecoding",
         path = res[[2 ;;, All, 2]];
         If[ model === "depModel",
             path[[All, All, ;; -2, 2]] = Map[diplotohaplo[[#]] &, path[[All, All, ;; -2, 2]], {2}];
         ];
         path = Map[CtStringFormat[#, "-"] &, path, {2}];
         path = StringJoin[Riffle[#, "||"]] & /@ path;
         path = Join[{{"Lines", "ViterbiPath"}}, Transpose[{sampleid, path}]];
         (*export*)
         summary = Join[{{key,"Genetic map of biallelic markers"}}, snpmap, 
             {{key, "Ln marginal likelihood"}}, logl2,
             {{key, "haplotypes in order"}}, haplotypes2,
             {{key, "diplotypes in order"}}, diplotypes2,
           {{key, "Viterbi path of "<>If[ model === "depModel",
                                          "haplotypes",
                                          "diplotypes"
                                      ]}}, path];
         Export[summaryFile, ExportString[summary, "CSV"], "Table"],
         "origPathSampling",
         path = res[[2 ;;, All, 2]];
         If[ model === "depModel",
             path[[All, All, All, ;; -2, 2]] = Map[diplotohaplo[[#]] &, path[[All, All, All, ;; -2, 2]], {3}];
         ];
         path = Map[CtStringFormat[#, "-"] &, path, {3}];
         path = Transpose[#] & /@ path;
         path = Map[StringJoin[Riffle[#, "|"]] &, path, {2}];
         path = Map[StringJoin[Riffle[#, "||"]] &, path, {1}];
         path = Join[{{"Lines", "SampledPaths"}}, Transpose[{sampleid, path}]];
         (*export*)
         summary = Join[{{key,"Genetic map of biallelic markers"}}, snpmap, 
           {{key, "Ln marginal likelihood"}}, logl2,
           {{key, "haplotypes in order"}}, haplotypes2,
           {{key, "diplotypes in order"}}, diplotypes2,
           {{key, "Sampled paths of "<>If[ model === "depModel",
                                           "haplotypes",
                                           "diplotypes"
                                       ]}}, path];
         Export[summaryFile, ExportString[summary, "CSV"], "Table"],
         _,
         Print["Wrong " <> resultFile <> "!"];
         Abort[]
         ]
    ]    

(*http://stackoverflow.com/questions/7525782/import-big-files-arrays-with-mathematica*)
readTable[file_String?FileExistsQ, format_String, chunkSize_: 1000] :=
    Module[ {stream, dataChunk, result, linkedList, add},
        SetAttributes[linkedList, HoldAllComplete];
        add[ll_, value_] :=
            linkedList[ll, value];
        stream = StringToStream[Import[file, "String"]];
        Internal`WithLocalSettings[Null,(*main code*)
            result = linkedList[];
            While[dataChunk =!= {}, 
             dataChunk = ImportString[StringJoin[Riffle[ReadList[stream, "String", chunkSize], "\n"]], format];
             result = add[result, dataChunk];
            ];
            result = Flatten[result, Infinity, linkedList],(*clean-up*)
         Close[stream]
         ];
        Join @@ result
    ]
    
getSummaryMR[summaryFile_String?FileExistsQ] :=
    Module[ {res, description,key = "MagicReconstruct-Summary"},
        res = readTable[summaryFile, "CSV"];
        res = Partition[Split[res, #1[[1]] != key && #2[[1]] != key &], 2];
        description = Flatten[#] & /@ res[[All, 1]];
        res = Join[{description}, res[[All, 2]]];
        res
    ]
     
End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicReconstruct`*"];

EndPackage[]


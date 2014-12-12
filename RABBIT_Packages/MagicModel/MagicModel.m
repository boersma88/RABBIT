(* Mathematica Package *)

(* Created by the Wolfram Workbench 30-Nov-2014 *)

BeginPackage["MagicModel`",{"MagicOrigin`","MagicOriginXY`"}]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["MagicModel`*"];
ClearAll @@ Names["MagicModel`*"];

magicPriorProcess::usage = "magicPriorProcess  "

magicSNPLikelihood::usage = "magicSNPLikelihood  "

Begin["`Private`"]
(* Implementation of the package *)

   
jointTransitionProb2 = 
  Function[{deltt, transitionRate, transitionRate2},
   Module[ {temp},
       temp = deltt transitionRate + deltt^2/2 transitionRate2;
       temp += IdentityMatrix[Length[temp]];
       temp
   ]
   ]

(*effective generation= {R^m,R^p} 
for maternally and paternally derived chromosomes.
R^m = R^p for autosomes*)
indepTransitionProb = Function[{deltt, nFounder, mapR},
  Module[ {mp,aa, aals},
      aals = Table[
        aa = Exp[-mp deltt] IdentityMatrix[nFounder];
        aa += (1 -IdentityMatrix[nFounder]) (1 - Exp[-mp deltt])/(nFounder -1), {mp, mapR}];
      KroneckerProduct @@ aals
  ]
  ]

(*mapR =(R^m+R^p)/2 for autsomes*)
depTransitionProb = Function[{deltt, nFounder, mapR},
   Module[ {ii, aa, bb},
       aa = Exp[-mapR deltt] IdentityMatrix[nFounder];
       aa += (1 - IdentityMatrix[nFounder]) (1-Exp[-mapR deltt])/(nFounder - 1);
       bb = ConstantArray[0, {nFounder^2, nFounder^2}];
       ii = (Range[nFounder] - 1) nFounder + Range[nFounder];
       bb[[ii, ii]] = aa;
       bb
   ]
   ]
   
(*mapR =R^m for the maternally derived X chromosome of a male*)
maleXTransitionProb = Function[{deltt, nFounder, mapR},
   Module[ {aa},
       aa = Exp[-mapR deltt] IdentityMatrix[nFounder];
       aa += (1 - IdentityMatrix[nFounder]) (1-Exp[-mapR deltt])/(nFounder - 1);
       aa
   ]
   ]


calMarkovProb[deltd_, model_, nFounder_,matingScheme_] :=
    Module[ {inbred, junction, mapR, startProb, rr, transitionRate, 
      transitionRate2, tranProb},
        (*{inbred, mapR, sumRho,junction}*)
        {inbred, mapR, junction} = N[origSummary[nFounder,matingScheme][[{1,2,4}]]];
        Switch[model,
         "jointModel",
         startProb = N[magicStationaryProb[nFounder, inbred[[-1]]]];
         rr = magicBasicRate[nFounder, inbred[[-1]], junction[[All, -1]]];
         transitionRate = magicRateMatrix[nFounder, rr];
         transitionRate2 = transitionRate.transitionRate;
         tranProb = Map[jointTransitionProb2[#, transitionRate, transitionRate2] &, deltd, {2}],
         "indepModel",
         startProb = Table[1./nFounder^2, {nFounder^2}];
         tranProb = Map[indepTransitionProb[#, nFounder, {mapR[[-1]],mapR[[-1]]}] &, deltd, {2}],
         "depModel",
         startProb = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1];
         startProb = startProb /. {{x_, x_} :> 1./nFounder, {_, _} -> 0};
         tranProb = Map[depTransitionProb[#, nFounder, mapR[[-1]]] &, deltd, {2}];];
        {startProb, tranProb}
    ]
       
calMarkovProbFemale[deltd_, model_, nFounder_, matingScheme_] :=
    Module[ {inbred, junction, mapR,startProb, rr, transitionRate, 
      transitionRate2, tranProb},
        (*{inbred, mapR, sumRho,junction}*)
        (*mapR: mapR[[1]] refers to R^m at all generations,and mapR[[2]] refers to R^p at all generations*)
        {inbred, mapR, junction} = N[origXYSummary[nFounder,matingScheme][[{1,2,4}]]];
        Switch[model,
          "jointModel",
          startProb = magicStationaryProbXY[nFounder, inbred[[-1]]];
          rr = magicBasicRateXY[nFounder, inbred[[-1]], junction[[All, -1]]];
          transitionRate = magicRateMatrixXY[nFounder, rr];
          transitionRate2 = transitionRate.transitionRate;
          tranProb = Map[jointTransitionProb2[#, transitionRate, transitionRate2] &, deltd,{2}],
          "indepModel",
          startProb = Table[1./nFounder^2, {nFounder^2}];
          tranProb = Map[indepTransitionProb[#, nFounder, mapR[[All,-1]]] &, deltd,{2}],
          "depModel",
          startProb = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1];
          startProb = startProb /. {{x_, x_} :> 1./nFounder, {_, _} -> 0};
          tranProb = Map[depTransitionProb[#, nFounder, Mean[mapR[[All,-1]]]] &, deltd,{2}];        
        ];
        {startProb, tranProb}
    ]  

(*only model like "depModel" is meaningfull for the maternally derived X chromsome of a male*)
calMarkovProbMale[deltd_,  nFounder_,matingScheme_] :=
    Module[ {mapR,startProb, tranProb},
        (*mapR: mapR[[1]] refers to R^m at all generations,and mapR[[2]] refers to R^p at all generations*)
        mapR = N[origXYSummary[nFounder,matingScheme][[2]]];
        startProb = Table[1./nFounder, {nFounder}];
        (*mapR[[1,-1]] take the R^m of the last generation*)
        tranProb = Map[maleXTransitionProb[#, nFounder, mapR[[1,-1]]] &, deltd,{2}];
        {startProb, tranProb}
    ]  

   
magicPriorProcess[geneticMap_,model_,nFounder_,matingScheme_] :=
    Module[ {map,chrs,deltd,posX,posA,startProbAutosome, tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale},
        startProbAutosome = tranProbAutosome = startProbFemale = tranProbFemale = startProbMale = tranProbMale = 0;
        map = SplitBy[Transpose[geneticMap[[2;;,2;;]]],First];
        (*change the genetic distance from centiMorgan into Morgan*)
        deltd = 0.01 Differences[#]&/@map[[All,All,2]];
        chrs = map[[All,1,1]];
        posX = Flatten[Position[chrs, "X"]];
        posA = Complement[Range[Length[chrs]], posX];
        If[ posA=!={},
            {startProbAutosome, tranProbAutosome} = calMarkovProb[deltd[[posA]],model, nFounder,matingScheme];
        ];
        If[ posX=!={},
            {startProbFemale, tranProbFemale} = calMarkovProbFemale[deltd[[posX]], model, nFounder,matingScheme];
            {startProbMale, tranProbMale} = calMarkovProbMale[deltd[[posX]], nFounder,matingScheme];
        ];
        {startProbAutosome, tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale}
    ]   



       
(*return res[[Y,IBD or not, D]] the probability of observed Y given the derived D and IBD status*)
likelihoodDiplo[epsF_, eps_] :=
    Module[ {likeli, likeliF,likeliFibd0, likeliFibd1, posterior},
     (*There are only two alleles 1 and 2, and missing alllels are denoted by N(Null).*)
     (*likeli, dimension 6 x4, refers to the probability of observed genotype given latent genotpe*)
     (*The 6 observed genotype is defined as {"NN" -> 1, "N1" -> 2, "1N" -> 2, "N2" -> 3, "2N" -> 3, "11" -> 4, "12" -> 5, "21" -> 5, "22" -> 6}     *)
     (*likeliFibd0 or likeliFibd1, dimension 9 x4, refers to the probability of derived genotype given latent genotype*)
     (*The 9 derived genotypes are defined as {"NN" -> 1, "N1" -> 2, "1N" -> 3, "N2" -> 4, "2N" -> 5, "11" -> 6, "12" -> 7, "21" -> 8, "22" -> 9}*)
     (*assuming that prior for latent genotype (1,1), (1,2), (2,1), and (22) are equally probable*)
        likeliFibd0 = {{1, 1, 1, 1}, 
                       {1 - epsF, epsF, 1 - epsF, epsF}, 
                       {1 - epsF, 1 - epsF, epsF, epsF}, 
                       {epsF, 1 - epsF, epsF, 1 - epsF}, 
                       {epsF, epsF, 1 - epsF, 1 - epsF}, 
                       {(1 - epsF)^2, (1 - epsF) epsF, (1 - epsF) epsF, epsF^2}, 
                       {(1 - epsF) epsF, (1 - epsF)^2, epsF^2, (1 - epsF) epsF}, 
                       {(1 - epsF) epsF, epsF^2, (1 -epsF)^2, (1 - epsF) epsF}, 
                       {epsF^2, (1 - epsF) epsF, (1 - epsF) epsF, (1 - epsF)^2}};
        likeliFibd1 = {{1, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
                       {0, 0, 0, 0}, {0, 0, 0, 0}, {1 - epsF, 0, 0, epsF}, 
                       {0, 0, 0, 0}, {0, 0, 0, 0}, {epsF, 0, 0, 1 - epsF}};
        likeli = {{1, 1, 1, 1}, {1 - eps, 1/2, 1/2, eps}, {eps, 1/2, 1/2, 1 - eps}, 
                  {(1 - eps)^2, (1 - eps) eps, (1 - eps) eps, eps^2}, 
                  {2 (1 - eps) eps, (1 - eps)^2 + eps^2, (1 - eps)^2 + eps^2, 2 (1 - eps) eps}, 
                  {eps^2, (1 - eps) eps, (1 - eps) eps, (1 - eps)^2}};
        Transpose[Simplify[Table[
           posterior = Transpose[Normalize[#, Total] & /@ likeliF];
           likeli.posterior, {likeliF, {likeliFibd0, likeliFibd1}}]]]
    ]

(*return res[[Y,D]] the probability of observed Y given the derived D *)
likelihoodHaplo[epsF_, eps_] :=
    Module[ {likeli, likeliF,posterior},
     (*There are only two alleles 1 and 2, and missing alllels are denoted by N(Null).*)
     (*likeli, dimension 3x2, refers to the probability of observed allele given latent true allele*)
     (*The 3 observed allele is defined as {"N" -> 1, "1" -> 2, "2" -> 3}*)
     (*likeliF, dimension 3x2,, refers to the probability of derived allele given latent true allele*)
     (*The 3 derived genotypes are defined as {"N" -> 1, "1" -> 2, "2" -> 3}*)
     (*assuming that prior for the two latent true alleles  1 and 2 are equally probable*)
        likeliF = {{1, 1},{1-epsF,epsF},{epsF,1-epsF}};
        likeli = {{1, 1},{1-eps,eps},{eps,1-eps}};
        posterior = Transpose[Normalize[#, Total] & /@ likeliF];
        likeli.posterior
    ]
        
magicSNPLikelihood[founderHaplo_, obsGeno_, epsF_, eps_, posA_, posX_, gender_] :=
    Module[ {likelidiplo,likelihaplo,posDiplo, posHaplo, nFounder, states, ibd12, 
      flatstates, derivedRule, obsRule, obsSeq, dataProb, hap, der, i, j},
      	likelidiplo = likelihoodDiplo[epsF, eps];
        likelihaplo = likelihoodHaplo[epsF, eps];
        {posDiplo, posHaplo} = 
         If[ gender == "Male",
             {posA, posX},
             {Join[posA, posX], {}}
         ];
        (*There are only two alleles 1 and 2,and missing alllels are denoted by N(Null).*)
        (*The second part of obsRule is for the X chromosome of a male*)
        obsRule = Join[{"NN" -> 1, "N1" -> 2, "1N" -> 2, "N2" -> 3, "2N" -> 3, 
           "11" -> 4, "12" -> 5, "21" -> 5, "22" -> 6}, {"N" -> 1, "1" -> 2,"2" -> 3}];
        obsSeq = Replace[obsGeno, obsRule, {2}];
        dataProb = Table[0, {Length[posDiplo] + Length[posHaplo]}];
        nFounder = Length[founderHaplo];
        (*caluclate dataProb for diplo*)
        If[ posDiplo=!={},
            states = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1];
            ibd12 = Boole[Equal @@ # & /@ states] + 1;
            flatstates = Flatten[states];
            derivedRule = {{"N", "N"} -> 1, {"N", "1"} -> 2, {"1", "N"} -> 
               3, {"N", "2"} -> 4, {"2", "N"} -> 5, {"1", "1"} -> 
               6, {"1", "2"} -> 7, {"2", "1"} -> 8, {"2", "2"} -> 9};
            dataProb[[posDiplo]] = Table[
              hap = founderHaplo[[All, i, j]];
              der = Replace[Partition[hap[[flatstates]], 2], derivedRule, {1}];
              Extract[likelidiplo[[obsSeq[[i, j]]]], Transpose[{ibd12, der}]]
              , {i, posDiplo}, {j, Dimensions[founderHaplo[[All, i]]][[2]]}];
        ];
        (*caluclate dataProb for haplo*)
        If[ posHaplo=!={},
            dataProb[[posHaplo]] = Table[
                der = founderHaplo[[All, i, j]] /. {"N" -> 1, "1" -> 2, "2" -> 3};
                likelihaplo[[obsSeq[[i, j]], der]], {i, posHaplo}, {j, Dimensions[founderHaplo[[All, i]]][[2]]}];
        ];
        dataProb
    ]

   

End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicModel`*"];

EndPackage[]


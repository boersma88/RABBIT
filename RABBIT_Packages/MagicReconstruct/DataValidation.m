(* Mathematica Package *)

BeginPackage["MagicReconstruct`DataValidation`"]

Unprotect @@ Names["MagicReconstruct`DataValidation`*"];
ClearAll @@ Names["MagicReconstruct`DataValidation`*"];

SNPValidation::usage = "SNPValidation  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 


checkgeneticmap[magicSNP_] :=
    Module[ {ls, ls2},
        ls = Transpose[magicSNP[[3 ;; 4, 2 ;;]]];
        ls = SplitBy[ls, First];
        ls2 = ls[[All, All, 2]];
        ls2 = OrderedQ[#] & /@ ls2;
        ls2 = Pick[ls[[All, 1, 1]], ls2, False];
        If[ ls2 == {},
            True,
            Print["The genetic distances must be strictly increasing on chromosomes ", ls2];
            False
            ]
    ]

checkfounderSNP[magicSNP_] :=
    Module[ {nfounder, ls},
        nfounder = magicSNP[[1, 2]];
        If[ ! (IntegerQ[nfounder] && nfounder >= 2),
            Print["the number of inbred founders must be an integer greater than or equal to 2!"]
        ];
        ls = Transpose[magicSNP[[5 ;; 4 + nfounder, 2 ;;]]];
        ls = ToString[#] & /@ Union[Flatten[ls]];
        ls = Complement[ls, {"1", "2", "N"}];
        If[ ls == {},
            True,
            Print["The alleles of founder haplotypes must be {1,2,N},where N denotes missing alelle.\nThe alleles ", ls, 
             " are not allowed! Make sure that the first founder haplotype is in row 5 and the number of founders is correct!"];
            False
        ]
    ]

checksampleSNP[magicSNP_] :=
    Module[ {nfounder, chrs, chr, posA, posX, geno, ls, pos, ls1, ls2, 
      res = True},
        nfounder = magicSNP[[1, 2]];
        chrs = Split[magicSNP[[3, 2 ;;]]][[All, 1]];
        chr = SplitBy[magicSNP[[3, 2 ;;]]][[All, 1]];
        If[ ! (Length[chr] == Length[Union[chr]]),
            Print["The SNP markers on the same chromsomes must be in order!"];
            res = False
        ];
        posX = Flatten[Position[chrs, "X"]];
        posA = Complement[Range[Length[chrs]], posX];
        geno = SplitBy[Transpose[Join[magicSNP[[{3}, 2 ;;]], magicSNP[[5 + nfounder ;;, 2 ;;]]]], First];
        ls = Transpose[Transpose[#] & /@ geno[[posA, All, 2 ;;]]];
        ls = Map[ToString, Union[Flatten[#]] & /@ ls, {2}];
        ls = Complement[#, {"NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"}] & /@ ls;
        pos = Flatten[Position[ls, _?(# =!= {} &), {1}, Heads -> False]];
        If[ pos =!= {},
            Print["The genotypes for sampled individuals on autosomals must be ", {"NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"}, ".\n", 
             "The following sampled individuals have illegal genotypes: \n", 
             Join[{{"Sampled individual", "Illegal genotypes"}}, 
               Transpose[{magicSNP[[5 + nfounder ;;, 1]][[pos]], ls[[pos]]}]] //TableForm];
            res = False;
        ];
        ls = Transpose[Transpose[#] & /@ geno[[posX, All, 2 ;;]]];
        ls = Map[ToString, Union[Flatten[#]] & /@ ls, {2}];
        ls1 = Complement[#, {"NN", "N1", "1N", "N2", "2N", "11", "12", "21",
               "22"}] === {} & /@ ls;
        ls2 = Complement[#, {"N", "1", "2"}] === {} & /@ ls;
        pos = Flatten[Position[MapThread[Or, {ls1, ls2}], False]];
        If[ pos =!= {},
            Print["The genotypes on female XX chromsomes must be ", {"NN", 
              "N1", "1N", "N2", "2N", "11", "12", "21", "22"}, ".\n", 
             "The haplotypes on male X chromsomes must be ", {"N", "1", "2"}, 
             ".\n", 
             "The following sampled individuals have illegal genotypes on X chromosomes: \n", 
             Join[{{"Sampled individual", "Illegal genotypes"}}, 
               Transpose[{magicSNP[[5 + nfounder ;;, 1]][[pos]], ls[[pos]]}]] //TableForm];
            res = False;
        ];
        res
    ]

SNPValidation[magicSNP_] :=
    If[ ! (checkgeneticmap[magicSNP] && checkfounderSNP[magicSNP] && 
        checksampleSNP[magicSNP]),
        Abort[]
    ]
     

End[] (* End Private Context *)

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicReconstruct`DataValidation`*"];

EndPackage[]
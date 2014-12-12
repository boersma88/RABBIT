(* Mathematica Package *)

BeginPackage["MaginOrigin`DOOrigin`", {"MagicOrigin`LPOrigin`"}]

inbredProbDO::usage = "inbredProbDO  "

mapExpansionDO::usage = "mapExpansionDO  "

juncDensityDO::usage = "juncDensityDO  "

sumJuncDensityDO::usage = "sumJuncDensityDO  "
(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

(*initial conditions refer to the F1 populations of DO*)
(*inialpha12 ~ beta_0(12)*)
(*inialpha123  ~ beta_0(123)*)    
(*iniR =R_0+alpha_0(12)*)
(*iniJK1122~0*)
(*iniJK1232 ~ K_0(1232)+alpha_0(123); 
*)

getinialpha12[nStrain_] :=1 - 1/nStrain  
    
getinialpha123[nStrain_] := (1 - 1/nStrain) (1 - 2/nStrain)

getiniR[founderAlpha12_,founderR_] :=founderR +founderAlpha12
    
getiniJK1122[] :=0
       
getiniJK1232[founderAlpha12_,founderR_,nStrain_] :=founderR (1-2/nStrain)+founderAlpha12 (1 - 2/nStrain)

inbredProbDO[founderAlpha12_,nStrain_,coalProb_,gCross_] :=
    Module[ {inialpha12,a12,t},
        inialpha12 = getinialpha12[nStrain];
        a12 = Table[Evaluate[lpnonIBDProb2[inialpha12,coalProb,t]],{t,gCross}];
        Join[{1-founderAlpha12,1-inialpha12},1 - a12]
    ]    

mapExpansionDO[founderAlpha12_,founderR_,nStrain_,coalProb_,gCross_] :=
    Module[ {inialpha12,iniR,Ra,t},
        inialpha12 = getinialpha12[nStrain];
        iniR = getiniR[founderAlpha12,founderR];
        Ra = Table[Evaluate[lpmapExpansion[inialpha12, iniR, coalProb,t]],{t,gCross}];
        Join[{founderR,iniR},Ra]
    ]       
    
junc1122DO[founderAlpha12_,founderR_,founderJ1122_, nStrain_,coalProb_, gCross_] :=
    Module[ {inialpha12,iniR, iniJK1122,j1122,t},
        inialpha12 = getinialpha12[nStrain];
        iniR = getiniR[founderAlpha12,founderR];
        iniJK1122 = getiniJK1122[];
        j1122 = Table[Evaluate[lpjunc1122[inialpha12, iniR, iniJK1122, coalProb, t]],{t,gCross}];
        Join[{founderJ1122,iniJK1122}, j1122]
    ]   
    
junc1232DO[founderAlpha12_,founderR_,founderJ1232_,nStrain_,coalProb_, gCross_] :=
    Module[ {inialpha12, inialpha123,iniJK1232,j1232,t},
        inialpha12 = getinialpha12[nStrain];
        inialpha123 = getinialpha123[nStrain]; 
        iniJK1232 = getiniJK1232[founderAlpha12,founderR,nStrain];
        j1232 = Table[Evaluate[lpjunc1232[inialpha123, iniJK1232,coalProb,t]],{t,gCross}];
        Join[{founderJ1232,iniJK1232},j1232]
    ]    

juncDensityDO[founderAlpha12_,founderR_,founderJ1122_,founderJ1232_,nStrain_,coalProb_, gCross_] :=
    Module[ {Ra,j1122,j1222,j1232},
        Ra = mapExpansionDO[founderAlpha12,founderR,nStrain,coalProb,gCross];
        j1122 = junc1122DO[founderAlpha12,founderR,founderJ1122, nStrain,coalProb, gCross];
        j1232 = junc1232DO[founderAlpha12,founderR,founderJ1232,nStrain,coalProb, gCross];
        j1222 = (Ra - j1232 - j1122)/2;
        {j1122,j1222,j1232}
    ]

sumJuncDensityDO[founderAlpha12_,founderR_,founderJ1122_, nStrain_,coalProb_, gCross_]  :=
    Module[ {Ra,j1122},
        Ra = mapExpansionDO[founderAlpha12,founderR,nStrain,coalProb,gCross];
        j1122 = junc1122DO[founderAlpha12,founderR,founderJ1122, nStrain,coalProb, gCross];
        2 Ra - j1122
    ]

End[] (* End Private Context *)

EndPackage[]
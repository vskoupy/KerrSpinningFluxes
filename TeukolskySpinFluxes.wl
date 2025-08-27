(* ::Package:: *)

(* ::Section::Closed:: *)
(*Begin package*)


BeginPackage["TeukolskySpinFluxes`",
  {"KerrGeodesics`",
  (*"Teukolsky`",*)
  "SpinningOrbit`"}
];


TeukolskySpinModeCircularCorrectionAnalytical::usage = "TeukolskySpinModeCircularCorrectionAnalytical[l,m,orbit,{angparNew,RCorrection}]"

TeukolskySpinModeEquatorialCorrectionAnalytical::usage = "TeukolskySpinModeEquatorialCorrectionAnalytical[l,m,n,orbit,{angparNew,RCorrection}]"
TeukolskySpinModeEquatorialCorrectionAnalyticalNew::usage = "TeukolskySpinModeEquatorialCorrectionAnalyticalNew[l,m,n,orbit,{angparNew,RCorrection}]"
TeukolskySpinModeEquatorialCorrectionAnalyticalNew2::usage = "TeukolskySpinModeEquatorialCorrectionAnalyticalNew2[l,m,n,orbit,{angparNew,RCorrection}]"

TeukolskySpinMode::usage = "TeukolskySpinMode[l,m,n,k,orbit] calculates Teukolsky amplitudes and fluxes from a spinning orbit";
TeukolskySpinModeFromCorrection::usage = "TeukolskySpinModeFromCorrection[l,m,n,k,orbitCorrection,s] calculates Teukolsky amplitudes and fluxes from a linear correction to an orbit and a value of the spin";
TeukolskySpinModeFromCorrectionNew::usage = "TeukolskySpinModeFromCorrectionNew[l,m,n,k,orbitCorrection,s] calculates Teukolsky amplitudes and fluxes from a linear correction to an orbit and a value of the spin";

TeukolskySpinModeAnalytical::usage = "TeukolskySpinModeFromCorrectionAalytical[l,m,n,k,orbit] calculates Teukolsky amplitudes and fluxes from analytical orbit";
TeukolskySpinModeAnalyticalNew::usage = "TeukolskySpinModeFromCorrectionAalyticalNew[l,m,n,k,orbit] calculates Teukolsky amplitudes and fluxes from analytical orbit"
TeukolskySpinModeAnalyticalNew2::usage = "TeukolskySpinModeFromCorrectionAalyticalNew2[l,m,n,k,orbit] calculates Teukolsky amplitudes and fluxes from analytical orbit"

TeukolskySpinModeCorrectionNum::usage = "TeukolskySpinModeCorrectionNum[l,m,n,k,orbitCorrection,\[Delta]\[Omega]] calculates linear correction to the fluxes from numerical derivatives of R ans S with step \[Delta]\[Omega]";
TeukolskySpinModeCorrection::usage = "TeukolskySpinModeCorrection[l,m,n,k,orbitCorrection,{angparNew,RCorrection}] calculates linear correction to the fluxes using given function for the R and S derivatives";
TeukolskySpinModeCorrectionAnalytical::usage = "TeukolskySpinModeCorrectionAnalytical[l,m,n,k,orbit,{angparNew,RCorrection}]"

TeukolskySpinModeSpherical::usage = "TeukolskySpinModeSpherical[l,m,k,orbitCorrection,s] calculates Teukolsky amplitudes and fluxes from a correction to a spherical orbit and a value of the spin";
TeukolskySpinModeSphericalCorrectionNum::usage = "TeukolskySpinModeSphericalCorrectionNum[l,m,k,orbitCorrection,\[Delta]\[Omega]] calculates linear correction to the fluxes for spherical orbit from numerical derivatives of R ans S with step \[Delta]\[Omega]";
TeukolskySpinModeSphericalCorrection::usage = "TeukolskySpinModeSphericalCorrection[l,m,k,orbitCorrection,orbitDerivatives,{angparNew,TeukolskySolverHS1spin}] calculates linear correction to the fluxes and their r and x derivatives for spherical orbit using given function for the R and S derivatives";
TeukolskySpinModeSphericalCorrectionNew::usage = "TeukolskySpinModeSphericalCorrectionNew[l,m,k,orbitCorrection,orbitDerivatives,{angparNew,TeukolskySolverHS1spin}] calculates linear correction to the fluxes and their r and x derivatives for spherical orbit using given function for the R and S derivatives";
TeukolskySpinModeSphericalCorrectionAnalytical::usage = "TeukolskySpinModeSphericalCorrectionAnalytical[l,m,k,orbitCorrection,orbitDerivatives,{angparNew,TeukolskySolverHS1spin}] calculates linear correction to the fluxes and their r and x derivatives for spherical orbit using given function for the R and S derivatives";

Trajectory::usage = "Trajectory[orbit] returns the orbital functions for orbit = KerrGeoOrbit[a,p,e,x]";


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Spinning fluxes*)


(* ::Subsection::Closed:: *)
(*Functions f_ab^(i), F_ab*)


fabi[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,a_]:={-2*\[Zeta]^2/\[CapitalDelta]^2*(L1L2S-2I*a/\[Zeta]*sin\[Theta]*L2S),
                                                       4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*((I*K/\[CapitalDelta]+1/\[Zeta]+1/\[Zeta]bar)*L2S-a*sin\[Theta]*K/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S),
                                                       4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*(L2S+I*a*sin\[Theta]*(1/\[Zeta]bar-1/\[Zeta])*S),
                                                       \[Zeta]^2/\[Zeta]bar^2*(I*(dK*\[CapitalDelta]-K*d\[CapitalDelta])-2/\[Zeta]*I*K*\[CapitalDelta]+K^2)/\[CapitalDelta]^2*S,
                                                       -2*\[Zeta]^2/\[Zeta]bar^2*(1/\[Zeta]+I*K/\[CapitalDelta])*S,
                                                       -\[Zeta]^2/\[Zeta]bar^2*S};
                                                       
dfabidS[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,K1_,dK1_,S_,L2S_,L1L2S_,S1_,L2S1_,L1L2S1_,a_]:={
      -2*\[Zeta]^2/\[CapitalDelta]^2*(L1L2S1-2I*a/\[Zeta]*sin\[Theta]*L2S1),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*((I*K/\[CapitalDelta]+1/\[Zeta]+1/\[Zeta]bar)*L2S1+(I*K1/\[CapitalDelta])*L2S-a*sin\[Theta]*K/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S1-a*sin\[Theta]*K1/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*(L2S1+I*a*sin\[Theta]*(1/\[Zeta]bar-1/\[Zeta])*S1),
      \[Zeta]^2/\[Zeta]bar^2*((I*(dK*\[CapitalDelta]-K*d\[CapitalDelta])-2/\[Zeta]*I*K*\[CapitalDelta]+K^2)/\[CapitalDelta]^2*S1+(I*(dK1*\[CapitalDelta]-K1*d\[CapitalDelta])-2/\[Zeta]*I*K1*\[CapitalDelta]+2*K*K1)/\[CapitalDelta]^2*S),
      -2*\[Zeta]^2/\[Zeta]bar^2*((1/\[Zeta]+I*K/\[CapitalDelta])*S1+(I*K1/\[CapitalDelta])*S),
      -\[Zeta]^2/\[Zeta]bar^2*S1};
      
dfabid\[Omega][\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dKdr_,dKd\[Omega]_,d2Kdrd\[Omega]_,S_,L2S_,L1L2S_,dSd\[Omega]_,dL2Sd\[Omega]_,dL1L2Sd\[Omega]_,a_]:={
      -2*\[Zeta]^2/\[CapitalDelta]^2*(dL1L2Sd\[Omega]-2I*a/\[Zeta]*sin\[Theta]*dL2Sd\[Omega]),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*((I*K/\[CapitalDelta]+1/\[Zeta]+1/\[Zeta]bar)*dL2Sd\[Omega]+(I*dKd\[Omega]/\[CapitalDelta])*L2S-a*sin\[Theta]*K/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*dSd\[Omega]-a*sin\[Theta]*dKd\[Omega]/\[CapitalDelta]*(1/\[Zeta]bar-1/\[Zeta])*S),
      4/Sqrt[2]*\[Zeta]^2/\[Zeta]bar/\[CapitalDelta]*(dL2Sd\[Omega]+I*a*sin\[Theta]*(1/\[Zeta]bar-1/\[Zeta])*dSd\[Omega]),
      \[Zeta]^2/\[Zeta]bar^2*((I*(dKdr*\[CapitalDelta]-K*d\[CapitalDelta])-2/\[Zeta]*I*K*\[CapitalDelta]+K^2)/\[CapitalDelta]^2*dSd\[Omega]+(I*(d2Kdrd\[Omega]*\[CapitalDelta]-dKd\[Omega]*d\[CapitalDelta])+2*(-1/\[Zeta]*I*\[CapitalDelta]+K)*dKd\[Omega])/\[CapitalDelta]^2*S),
      -2*\[Zeta]^2/\[Zeta]bar^2*((1/\[Zeta]+I*K/\[CapitalDelta])*dSd\[Omega]+(I*dKd\[Omega]/\[CapitalDelta])*S),
      -\[Zeta]^2/\[Zeta]bar^2*dSd\[Omega]};

dfabidr[\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,a_,\[Omega]_]:={
      1/\[CapitalDelta]^3 (\[CapitalDelta] (4 I a L2S sin\[Theta]-4 L1L2S \[Zeta])+4 \[Zeta] (-2 I a L2S sin\[Theta]+L1L2S \[Zeta]) d\[CapitalDelta]),
      1/(\[Zeta]bar^3 \[CapitalDelta]^3) 2 Sqrt[2] (a S (\[Zeta]-\[Zeta]bar) sin\[Theta] (2 d\[CapitalDelta] \[Zeta] \[Zeta]bar K-(dK \[Zeta] \[Zeta]bar+(-2 \[Zeta]+\[Zeta]bar) K) \[CapitalDelta])+L2S (-I \[Zeta] \[Zeta]bar K (2 d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-2 \[Zeta]bar) \[CapitalDelta])+\[CapitalDelta] (\[Zeta] \[Zeta]bar (I dK \[Zeta] \[Zeta]bar-d\[CapitalDelta] (\[Zeta]+\[Zeta]bar))-(\[Zeta]-\[Zeta]bar) (2 \[Zeta]+\[Zeta]bar) \[CapitalDelta]))),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) 2 Sqrt[2] (-L2S \[Zeta] \[Zeta]bar (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-2 \[Zeta]bar) \[CapitalDelta])-I a S (\[Zeta]-\[Zeta]bar) sin\[Theta] (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(2 \[Zeta]-\[Zeta]bar) \[CapitalDelta])),
      -(1/(\[Zeta]bar^3 \[CapitalDelta]^3))2 I S (-I \[Zeta] K^2 (d\[CapitalDelta] \[Zeta] \[Zeta]bar+(\[Zeta]-\[Zeta]bar) \[CapitalDelta])+\[Zeta]^2 \[CapitalDelta] (dK d\[CapitalDelta] \[Zeta]bar+(dK-\[Zeta]bar \[Omega]) \[CapitalDelta])+K (-d\[CapitalDelta]^2 \[Zeta]^2 \[Zeta]bar+\[Zeta]^2 (-d\[CapitalDelta]+\[Zeta]bar+I dK \[Zeta]bar) \[CapitalDelta]+(-2 \[Zeta]+\[Zeta]bar) \[CapitalDelta]^2)),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) 2 S (I \[Zeta] K (d\[CapitalDelta] \[Zeta] \[Zeta]bar+2 (\[Zeta]-\[Zeta]bar) \[CapitalDelta])+\[CapitalDelta] (-I dK \[Zeta]^2 \[Zeta]bar+(2 \[Zeta]-\[Zeta]bar) \[CapitalDelta])),
      (2 S \[Zeta] (\[Zeta]-\[Zeta]bar))/\[Zeta]bar^3};

dfabid\[Theta][\[Zeta]_,\[Zeta]bar_,sin\[Theta]_,\[CapitalDelta]_,d\[CapitalDelta]_,K_,dK_,S_,L2S_,L1L2S_,dSd\[Theta]_,dL2Sd\[Theta]_,dL1L2Sd\[Theta]_,a_]:={
      (-4 a^2 L2S sin\[Theta]^2+4 I a (dL2Sd\[Theta]-L1L2S) sin\[Theta] \[Zeta]+2 \[Zeta] (-((dL1L2Sd\[Theta]+L2S) \[Zeta])+L2S \[Zeta]bar))/\[CapitalDelta]^2,
      Sqrt[2] (-2 I a^2 K S sin\[Theta]^2 (2 \[Zeta]-\[Zeta]bar) (\[Zeta]+\[Zeta]bar)+\[Zeta] \[Zeta]bar (-I K S (\[Zeta]-\[Zeta]bar)^2+2 I dL2Sd\[Theta] K \[Zeta] \[Zeta]bar+2 dL2Sd\[Theta] \[CapitalDelta] (\[Zeta]+\[Zeta]bar))+2 a sin\[Theta] (dSd\[Theta] K \[Zeta] \[Zeta]bar (-\[Zeta]+\[Zeta]bar)+I L2S \[CapitalDelta] (\[Zeta]+\[Zeta]bar) (2 \[Zeta]+\[Zeta]bar)-K L2S \[Zeta] \[Zeta]bar (\[Zeta]+2 \[Zeta]bar)))/(\[CapitalDelta]^2 \[Zeta]bar^3),
      Sqrt[2] (-2 a^2 S sin\[Theta]^2 (2 \[Zeta]-\[Zeta]bar) (\[Zeta]+\[Zeta]bar)-\[Zeta] \[Zeta]bar (S (\[Zeta]-\[Zeta]bar)^2-2 dL2Sd\[Theta] \[Zeta] \[Zeta]bar)+2 I a sin\[Theta] \[Zeta] \[Zeta]bar (dSd\[Theta] (\[Zeta]-\[Zeta]bar)+L2S (\[Zeta]+2 \[Zeta]bar)))/(\[CapitalDelta] \[Zeta]bar^3),
      1/(\[Zeta]bar^3 \[CapitalDelta]^2) (\[Zeta] (dSd\[Theta] \[Zeta] \[Zeta]bar+2 I a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) K^2+dK \[Zeta] (I dSd\[Theta] \[Zeta] \[Zeta]bar-2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) \[CapitalDelta]+K (d\[CapitalDelta] \[Zeta] (-I dSd\[Theta] \[Zeta] \[Zeta]bar+2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar))+2 (-I dSd\[Theta] \[Zeta] \[Zeta]bar+a S sin\[Theta] (2 \[Zeta]+\[Zeta]bar)) \[CapitalDelta])),
      -(1/(\[Zeta]bar^3 \[CapitalDelta]))2 (\[Zeta] (I dSd\[Theta] \[Zeta] \[Zeta]bar-2 a S sin\[Theta] (\[Zeta]+\[Zeta]bar)) K+(dSd\[Theta] \[Zeta] \[Zeta]bar+I a S sin\[Theta] (2 \[Zeta]+\[Zeta]bar)) \[CapitalDelta]),
      -((\[Zeta] (dSd\[Theta] \[Zeta] \[Zeta]bar+2 I a S sin\[Theta] (\[Zeta]+\[Zeta]bar)))/\[Zeta]bar^3)};

Fab[\[Zeta]_,\[Zeta]bar_,a_,sin\[Theta]p_,RIn_,DRIn_,DDRIn_,S_,L2S_,L1L2S_] := {-\[Zeta]/(2*\[Zeta]bar)*(L1L2S - 2*I*a*sin\[Theta]p*L2S/\[Zeta])*RIn,
                                                             (sin\[Theta]p \[Zeta])/\[Zeta]bar*L2S*DRIn - (sin\[Theta]p (\[Zeta]+\[Zeta]bar))/\[Zeta]bar^2*L2S*RIn + (I a sin\[Theta]p^2 (\[Zeta]-\[Zeta]bar))/\[Zeta]bar^2*S*DRIn,
                                                             -\[Zeta]/(2 \[Zeta]bar)*S*sin\[Theta]p^2*(DDRIn - 2/\[Zeta]*DRIn),
                                                             1/2 I DRIn (a L2S sin\[Theta]p+I L1L2S \[Zeta]),
                                                             -(1/2) L2S sin\[Theta]p (DRIn-DDRIn \[Zeta])};

FabEq[r_,a_,RIn_,DRIn_,DDRIn_,S_,L2S_,L1L2S_] := {-1/2*(L1L2S - 2*I*a*L2S/r)*RIn,
                                                             L2S*DRIn - 2/r*L2S*RIn,
                                                             -1/2*S*(DDRIn - 2/r*DRIn),
                                                             1/2 I DRIn (a L2S+I L1L2S r),
                                                             -(1/2) L2S (DRIn-DDRIn r)};

dFabd\[Omega][\[Zeta]_,\[Zeta]bar_,a_,sin\[Theta]p_,RIn_,DRIn_,DDRIn_,dRInd\[Omega]_,dDRInd\[Omega]_,dDDRInd\[Omega]_,S_,L2S_,L1L2S_,dSd\[Omega]_,dL2Sd\[Omega]_,dL1L2Sd\[Omega]_] := {
    -\[Zeta]/(2*\[Zeta]bar)*((dL1L2Sd\[Omega] - 2*I*a*sin\[Theta]p*dL2Sd\[Omega]/\[Zeta])*RIn + (L1L2S - 2*I*a*sin\[Theta]p*L2S/\[Zeta])*dRInd\[Omega]),
    ((sin\[Theta]p \[Zeta])/\[Zeta]bar*dL2Sd\[Omega]*DRIn + (sin\[Theta]p \[Zeta])/\[Zeta]bar*L2S*dDRInd\[Omega] - (sin\[Theta]p (\[Zeta]+\[Zeta]bar))/\[Zeta]bar^2*dL2Sd\[Omega]*RIn - (sin\[Theta]p (\[Zeta]+\[Zeta]bar))/\[Zeta]bar^2*L2S*dRInd\[Omega] 
      + (I a sin\[Theta]p^2 (\[Zeta]-\[Zeta]bar))/\[Zeta]bar^2*dSd\[Omega]*DRIn + (I a sin\[Theta]p^2 (\[Zeta]-\[Zeta]bar))/\[Zeta]bar^2*S*dDRInd\[Omega]),
    -\[Zeta]/(2 \[Zeta]bar)*sin\[Theta]p^2*(dSd\[Omega]*(DDRIn - 2/\[Zeta]*DRIn) + S*(dDDRInd\[Omega] - 2/\[Zeta]*dDRInd\[Omega]))};

dFabd\[Omega]Eq[r_,a_,RIn_,DRIn_,DDRIn_,dRInd\[Omega]_,dDRInd\[Omega]_,dDDRInd\[Omega]_,S_,L2S_,L1L2S_,dSd\[Omega]_,dL2Sd\[Omega]_,dL1L2Sd\[Omega]_] := {
    (-r (dRInd\[Omega] L1L2S+dL1L2Sd\[Omega] RIn)+2 I a (dRInd\[Omega] L2S+dL2Sd\[Omega] RIn))/(2 r),
    (-2 dRInd\[Omega] L2S+dL2Sd\[Omega] DRIn r+dDRInd\[Omega] L2S r-2 dL2Sd\[Omega] RIn)/r,
    -((-2 DRIn dSd\[Omega]+DDRIn dSd\[Omega] r-2 dDRInd\[Omega] S+dDDRInd\[Omega] r S)/(2 r))}


(* ::Subsection::Closed:: *)
(*Circular*)


Options[TeukolskySpinModeCircularCorrectionAnalytical] = {WorkingPrecision->30};
TeukolskySpinModeCircularCorrectionAnalytical[l_?IntegerQ,m_?IntegerQ,orbit_,{angparNew_,RCorrection_},OptionsPattern[]]:=Module[{prec,a,p,e,x,
    En0,Lz0,K0,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t1,\[Omega],\[Omega]1,SWSH,dSWSHd\[Omega],R,Rp,\[Lambda],d\[Lambda]d\[Omega],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumPlus1,sumMinus0,sumMinus1,
    rp,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],V,dVd\[Omega],
    RIn,dRIndr,d2RIndr2,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,dDRInd\[Omega],d2DRIndrd\[Omega],dDDRInd\[Omega],dDRUpd\[Omega],d2DRUpdrd\[Omega],dDDRUpd\[Omega],
    \[Theta]p,K\[Theta],dK\[Theta]d\[Omega],S,dSd\[Theta],d2Sd\[Theta]2,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],
    L2S,dL2Sd\[Theta],L1L2S,dL2Sd\[Omega],d2L2Sd\[Theta]d\[Omega],dL1L2Sd\[Omega],
    vn,vmb,
    FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega],
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega],
    CPlus0,CPlus1,CMinus0,CMinus1,
    \[ScriptCapitalF]En\[ScriptCapitalI],\[ScriptCapitalF]En\[ScriptCapitalH],\[ScriptCapitalF]Lz\[ScriptCapitalI],\[ScriptCapitalF]Lz\[ScriptCapitalH],\[ScriptCapitalF]En\[ScriptCapitalI]1,\[ScriptCapitalF]En\[ScriptCapitalH]1,\[ScriptCapitalF]Lz\[ScriptCapitalI]1,\[ScriptCapitalF]Lz\[ScriptCapitalH]1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  prec = OptionValue[WorkingPrecision];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["Inclination"];
  If[x!=1,Return[$Failed]];
  En0 = orbit["Energy"]; (* Shifted constants of motion *)
  Lz0 = orbit["AngularMomentum"];
  K0 = orbit["CarterConstant"] + (Lz0-a*En0)^2;
  {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[orbit["Frequencies"]];(* Mino frequencies *)
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t; (* BL frequencies *)
  \[CapitalUpsilon]\[Tau] = orbit["ProperTimeFrequency"];
  \[CapitalUpsilon]t1 = -3*\[CapitalUpsilon]\[Tau]/(2*Sqrt[K0]);
  \[Omega] = m*\[CapitalOmega]\[Phi]; (* Frequency of mode *)
  \[Omega]1 = 3*\[CapitalUpsilon]\[Tau]*\[Omega]/(2*Sqrt[K0]*\[CapitalUpsilon]t);(* Linear part of the frequency *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,
                                     SetPrecision[a, prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega], prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                         "precODE" -> prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  (*Print["Calculated R: "<>ToString[AbsoluteTiming[*)
    R = RCorrection[-2,l,m,SetPrecision[a,prec+5],
                         SetPrecision[\[Omega],prec+5],1,
                         SetPrecision[\[Lambda],prec+5],
                         SetPrecision[d\[Lambda]d\[Omega],prec+5],e,p,"precODE"->prec];
  (*][[1]]]];*)
  \[Lambda]1 = d\[Lambda]d\[Omega]*\[Omega]1;
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  (* Functions of \[Theta] for equatorial orbits *)
  \[Theta]p = N[Pi/2,prec];
  K\[Theta] = a*\[Omega] - m;
  dK\[Theta]d\[Omega] = a;
  (* Spin-weighted spheroidal harmonics and their \[Theta]-derivatives *)
  {S,dSd\[Theta]}=SWSH[\[Theta]p]; 
  d2Sd\[Theta]2 = - (-a^2*\[Omega]^2 - m^2 - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
  (* \[Omega]-derivatives of S *)
  {dSd\[Omega],d2Sd\[Theta]d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta])  *)
  d3Sd\[Theta]2d\[Omega] = - (-a^2*\[Omega]^2 - m^2 - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega] + 2*m*a + d\[Lambda]d\[Omega])*S;
  (* Operators L and their \[Omega]-derivatives *)
  L2S = dSd\[Theta] + K\[Theta]*S;
  dL2Sd\[Theta] = d2Sd\[Theta]2 + K\[Theta]*dSd\[Theta] -2*S;
  L1L2S = dL2Sd\[Theta] + K\[Theta]*L2S; 
  dL2Sd\[Omega] = d2Sd\[Theta]d\[Omega] + K\[Theta]*dSd\[Omega] + dK\[Theta]d\[Omega]*S;
  d2L2Sd\[Theta]d\[Omega] = d3Sd\[Theta]2d\[Omega] + K\[Theta]*d2Sd\[Theta]d\[Omega] + dK\[Theta]d\[Omega]*dSd\[Theta] -2*dSd\[Omega];
  dL1L2Sd\[Omega] = d2L2Sd\[Theta]d\[Omega] + K\[Theta]*dL2Sd\[Omega] + dK\[Theta]d\[Omega]*L2S;
  rp = p;
  Rp = R[rp];(* Evaluate the radial functions and their derivatives *)
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKdr = 2*rp*\[Omega];
  dKd\[Omega]  = (rp^2+a^2);
  d2Kdrd\[Omega] = 2*rp;
  d\[CapitalDelta] = 2*(rp-1);
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  (* In solutions of the radial equation and their r and \[Omega] derivatives *)
  RIn    = Rp["R0"]["In"]["R"]; 
  dRIndr   = Rp["R0"]["In"]["dR"];
  d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  dRInd\[Omega]   = Rp["R1"]["In"]["R"]; (* \[Omega]-derivatives of the radial function *)
  d2RIndrd\[Omega]  = Rp["R1"]["In"]["dR"];
  d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRIn = dRIndr - I*K/\[CapitalDelta]*RIn;
  dDRInd\[Omega] = d2RIndrd\[Omega] - I*K/\[CapitalDelta]*dRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RIn;
  dDRIndr = d2RIndr2 - I*K/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  d2DRIndrd\[Omega] = d3RIndr2d\[Omega] - I*K/\[CapitalDelta]*d2RIndrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRInd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
  dDDRInd\[Omega] = d2DRIndrd\[Omega] - I*K/\[CapitalDelta]*dDRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRIn;
  (* Up solutions of the radial equation and their r and \[Omega] derivatives *)
  RUp    = Rp["R0"]["Up"]["R"]; 
  dRUpdr   = Rp["R0"]["Up"]["dR"];
  d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  dRUpd\[Omega]   = Rp["R1"]["Up"]["R"]; (* \[Omega]-derivative of the radial function *)
  d2RUpdrd\[Omega]  = Rp["R1"]["Up"]["dR"];
  d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRUp = dRUpdr - I*K/\[CapitalDelta]*RUp;
  dDRUpd\[Omega] = d2RUpdrd\[Omega] - I*K/\[CapitalDelta]*dRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RUp;
  dDRUpdr = d2RUpdr2 - I*K/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  d2DRUpdrd\[Omega] = d3RUpdr2d\[Omega] - I*K/\[CapitalDelta]*d2RUpdrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRUpd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
  dDDRUpd\[Omega] = d2DRUpdrd\[Omega] - I*K/\[CapitalDelta]*dDRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRUp;
  (* components of the velocity *)
  vn = -((rp^2+a^2)*En0 - a*Lz0)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
  vmb = (-I*(a*En0 - Lz0));
  (* In solution *)
  {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = FabEq[rp,a,RIn,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
  {dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RIn,DRIn,DDRIn,dRInd\[Omega],dDRInd\[Omega],dDDRInd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
  (* Up solution *)
  {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = FabEq[rp,a,RUp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
  {dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RUp,DRUp,DDRUp,dRUpd\[Omega],dDRUpd\[Omega],dDDRUpd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
  sumPlus0 = (vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn);
  sumPlus1 = (
      (vn*GnIn + vmb*GmbIn)/Sqrt[K0] +
      (vn*vn*dFnnInd\[Omega] + vn*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
   );
  sumMinus0 = (vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp);
  sumMinus1 = (
      (vn*GnUp + vmb*GmbUp)/Sqrt[K0] +
      (vn*vn*dFnnUpd\[Omega] + vn*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
   );
  W = (RIn*dRUpdr - dRIndr*RUp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  W1 = (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/(rp^2-2*rp+a^2)*\[Omega]1; (* derivative of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalUpsilon]t*W); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalUpsilon]t*W);
  CPlus1  = 2*Pi*sumPlus1/(\[CapitalUpsilon]t*W) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CPlus0; (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*sumMinus1/(\[CapitalUpsilon]t*W) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CMinus0;
  \[ScriptCapitalF]En\[ScriptCapitalI] = Abs[CPlus0]^2/(4Pi*\[Omega]^2); (* Fluxes and their linear parts *)
  \[ScriptCapitalF]En\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI] = Abs[CPlus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]En\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]En\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  <|
    "l"->l, "m"->m, "k"->0, "n"->0,
    "\[Omega]"->\[Omega], "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"-><|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"-><|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha], "\[Alpha]Correction"->\[Alpha]1,
    "S"->S,
    "Fluxes"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]|>
    |>,
    "FluxesCorrection"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]1|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]1|>
    |>
  |>
]


(*Options[TeukolskySpinModeCircularCorrection] = {WorkingPrecision->30};

TeukolskySpinModeCircularCorrection[l_,m_,orbitCorrection_,orbitDerivatives_,{angparNew_,TeukolskySolverHS1spin_},OptionsPattern[]]:=Module[{
    h1,h2,a,p,e,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,dEndr,dLzdr,dEndx,dLzdx,\[CapitalOmega]\[Phi],\[CapitalOmega]\[Phi]1,d\[CapitalOmega]\[Phi]dr,d\[CapitalOmega]\[Phi]dx,correction,derivatives,z,\[CapitalGamma],\[CapitalGamma]1,d\[CapitalGamma]dr,d\[CapitalGamma]dx,
    \[Omega],\[Omega]prec,\[Omega]1,d\[Omega]dr,d\[Omega]dx,\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega],R,\[ScriptCapitalC]2,d\[ScriptCapitalC]2d\[Omega],rplus,P,\[Epsilon],\[Alpha],d\[Alpha]d\[Omega],sumPlus0,sumMinus0,sumPlus1,sumMinus1,dsumPlusdr,dsumMinusdr,dsumPlusdx,dsumMinusdx,
    rp,\[CapitalDelta],d\[CapitalDelta],K,dKd\[Omega],dKdr,d2Kdrd\[Omega],V,dVd\[Omega],dVdr,
    RIn,dRIndr,d2RIndr2,d3RIndr3,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,d3RUpdr3,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    \[Theta]2,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],L2S,dL2Sd\[Theta],dL2Sd\[Omega],L1L2S,dL1L2Sd\[Theta],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0d\[Omega],dfnmb0d\[Omega],dfnmb1d\[Omega],dfmbmb0d\[Omega],dfmbmb1d\[Omega],dfmbmb2d\[Omega],dfnn0dr,
    dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],
    ul,un,um,umb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,
    St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,
    dzpdr,d\[Theta]pdr,dzpdx,d\[Theta]pdx,dUzpdr,dUzpdx,d\[CapitalSigma]dr,d\[CapitalSigma]dx,dexpdr,dexpdx,un1,dundr,dundx,umb1,dumbdr,dumbdx,
    Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,W,dWd\[Omega],
    CPlus0,CMinus0,CPlus1,CMinus1,dCPlusdr,dCMinusdr,dCPlusdx,dCMinusdx},
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_] := 1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+
               a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  Print["Calculating l = "<>ToString[l]<>", m = "<>ToString[m]<>", k = "<>ToString[k]<>" mode"];
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  If[e!=0,Return[$Failed]];
  \[ScriptCapitalI] = orbitCorrection["\[ScriptCapitalI]"];
  En = orbitCorrection["Ehat"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzhat"];
  Kc = orbitCorrection["Khat"];
  En1 = orbitCorrection["ES"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["LS"];
  dEndr = orbitDerivatives["dEndr"]; (* Derivatives of the constants of motion *)
  dLzdr = orbitDerivatives["dLzdr"];
  dEndx = orbitDerivatives["dEndx"]; 
  dLzdx = orbitDerivatives["dLzdx"];
  {\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Coordinate frequencies *)
  {\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the coordinate frequencies *)
  {d\[CapitalOmega]\[Phi]dr} = orbitDerivatives["dBLFrequenciesdr"]; (* Derivatives of the coordinate frequencies *)
  {d\[CapitalOmega]\[Phi]dx} = orbitDerivatives["dBLFrequenciesdx"];
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  derivatives = orbitDerivatives["OrbitDerivatives"]; (* function containing derivatives of the trajectory *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction and derivatives *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]]; 
  d\[CapitalGamma]dr = orbitDerivatives["dMinoFrequenciesdr"][[1]]; 
  d\[CapitalGamma]dx = orbitDerivatives["dMinoFrequenciesdx"][[1]]; 
  \[Omega] = m*\[CapitalOmega]\[Phi]; (* Geodesic frequency and the linear correction and derivatives *)
  \[Omega]prec = m*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[Cos[\[ScriptCapitalI]],OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"]; (* Frequency with higher precision *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1;
  d\[Omega]dr = m*d\[CapitalOmega]\[Phi]dr;
  d\[Omega]dx = m*d\[CapitalOmega]\[Phi]dx;
  (*Print["Calculating angular functions"<>ToString@AbsoluteTiming[{\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1];][[1]]];*)(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1,"precODE"->OptionValue[WorkingPrecision]];(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  Print["Calculating radial functions"<>ToString@AbsoluteTiming[R = TeukolskySolverHS1spin[p,-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1,\[Lambda],d\[Lambda]d\[Omega],"precODE"->OptionValue[WorkingPrecision]]][[1]]];
  (*R = TeukolskySolverHS1spin[p,-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1,\[Lambda],d\[Lambda]d\[Omega]];*)
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  d\[ScriptCapitalC]2d\[Omega] = 4 \[Lambda]^3 d\[Lambda]d\[Omega]+4 \[Lambda]^2 (3 d\[Lambda]d\[Omega]+10 a (m-2 a \[Omega]))+8 \[Lambda] (d\[Lambda]d\[Omega] (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega])) + 
          48 \[Omega] (a m d\[Lambda]d\[Omega]+6-18 a^3 m \[Omega]+12 a^4 \[Omega]^2+a^2 (d\[Lambda]d\[Omega] \[Omega]+6 m^2));  (* \[Omega]-derivative of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  d\[Alpha]d\[Omega] = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*d\[ScriptCapitalC]2d\[Omega] + 256*(2*rplus)^5*((P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2)/\[ScriptCapitalC]2; (* \[Omega]-derivative of the constant for horizon fluxes *)
  rp = p;
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKd\[Omega]  = (rp^2+a^2);
  d\[CapitalDelta] = 2*(rp-1);
  dKdr = 2*rp*\[Omega];
  d2Kdrd\[Omega] = 2*rp;
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  dVdr = -((2K*dKdr + 4I*K + 4I*(rp - 1)*dKdr)*\[CapitalDelta] - (K^2 + 4I*(rp - 1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2 + 8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
  {{{RIn,dRIndr},{dRInd\[Omega],d2RIndrd\[Omega]}},{{RUp,dRUpdr},{dRUpd\[Omega],d2RUpdrd\[Omega]}}} = R; (* Radial function and its derivatives *)
  d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from the Teukolsky equation *)
  d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  
  d3RIndr3 = 1/\[CapitalDelta] (dVdr*RIn + (V + 2)*dRIndr);  (* third derivative of radial function from the Teukolsky equation *)
  d3RUpdr3 = 1/\[CapitalDelta] (dVdr*RUp + (V + 2)*dRUpdr);
  d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
  d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
    {S,dSd\[Theta],d2Sd\[Theta]2}=SWSH[N[Pi/2,OptionValue[WorkingPrecision]]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
    {dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega]}=dSWSHd\[Omega][N[Pi/2,OptionValue[WorkingPrecision]]];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
    d3Sd\[Theta]3 = -(-(a*\[Omega])^2-(m)^2-2+2*m*a*\[Omega]+\[Lambda]-1)*dSd\[Theta]; (*third derivative from the TE *)
    L2S    = dSd\[Theta]-(S (m+a (-1) \[Omega])); (* Operators acting on S(\[Theta]) *)
    L1L2S  = d2Sd\[Theta]2+(dSd\[Theta] (-2 m-2 a (-1) \[Omega]))+S (-2-(m (m))/(-1)-2 a (m) \[Omega]-a^2 (-1) \[Omega]^2);  
    dL2Sd\[Omega] = (d2Sd\[Theta]d\[Omega]-(dSd\[Omega] (m+a (-1) \[Omega])+S (a (-1)))); (* \[Theta] and \[Omega]-derivatives of the operator *)
    dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1) (dSd\[Theta] (m+a (-1) \[Omega])+S (2));
    dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1) dSd\[Theta] (5-m^2-2 a (m) (-1) \[Omega]-a^2 \[Omega]^2)+
               (d2Sd\[Theta]2 (-1) (2 m+2 a (-1) \[Omega])+2 S (m+a \[Omega] (-2)));
    dL1L2Sd\[Omega] = (d3Sd\[Theta]2d\[Omega]+(d2Sd\[Theta]d\[Omega] (-2 m-2 a (-1) \[Omega])+dSd\[Theta] (-2 a (-1))) + 
               dSd\[Omega] (-2-(m (m))/(-1)-2 a (m) \[Omega]-a^2 (-1) \[Omega]^2) + S (-2 a (m) - a^2 (-1) 2*\[Omega]));  
    \[Zeta] = rp;
    \[Zeta]bar = rp;
    \[CapitalSigma] = rp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,a,\[Omega]];
    {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
    {dfnn0d\[Omega],dfnmb0d\[Omega],dfnmb1d\[Omega],dfmbmb0d\[Omega],dfmbmb1d\[Omega],dfmbmb2d\[Omega]} = dfabid\[Omega][\[Zeta],\[Zeta]bar,1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega],a];
    un = -((rp^2+a^2)*En - a*Lz)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
    ul = -((rp^2+a^2)*En - a*Lz)/(\[CapitalDelta]);
    um = (I*(a*En - Lz))/(-Sqrt[2]*\[Zeta]bar);
    umb = Conjugate[um];
    (* Dipole term *)
    Sln  = (-((rp (Kc))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
    Snm  = (\[Zeta]/Sqrt[Kc])*um*un;
    Snmb = Conjugate[Snm];
    Slmb = (-(\[Zeta]/Sqrt[Kc]))*ul*umb;
    Smmb = 0;
    Cmnn   = un^2;
    Cmnmb  = un*umb;
    Cmmbmb = umb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = 0;
    pi = -((I a)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-I a (-2))/(2 \[Zeta]^2 Sqrt[2]));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*un-Sln*2*Re[gamma]*un-2*Re[Snmb*((Conjugate[alpha]+beta)*un-mu*um)]);
    Cdmbmb= (Scd\[Gamma]mbdc*umb-Snmb*(-pi*ul)-Slmb*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)+Smmb*(-(-alpha+Conjugate[beta])*umb));
    Cdnmb = (Scd\[Gamma]ndc*umb+Scd\[Gamma]mbdc*un-Sln*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)-Snmb*(Conjugate[rho]*un-mu*ul-(Conjugate[alpha]-beta)*umb)
      -Snm*(-(-alpha+Conjugate[beta])*umb)-Snmb*(-Conjugate[pi]*umb-pi*um)-Slmb*(2*Re[gamma]*un)+Smmb*((alpha+Conjugate[beta])*un-Conjugate[mu]*umb))/2;
    St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]-m)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]-m)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
    S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
    (* Derivatives of the trajectory *)
    dzpdr = derivatives["dzdr"]; (* Derivatives of the coordinates and four-velocity *)
    dzpdx = derivatives["dzdx"];
    d\[Theta]pdr = -dzpdr;
    d\[Theta]pdx = -dzpdx;
    dUzpdr = derivatives["dUzdr"];
    dUzpdx = derivatives["dUzdx"];
    d\[CapitalSigma]dr = 2*rp;
    d\[CapitalSigma]dx = 0;
    dexpdr = I*(\[Omega]*derivatives["d\[CapitalDelta]tdr"]-m*derivatives["d\[CapitalDelta]\[Phi]dr"]);
    un1  = -( ((rp^2+a^2)*(En1-h1[rp,0]) - a*(Lz1+h2[rp,0]) ))/(2*\[CapitalSigma]);(* Linear parts and derivatives of the four-velocity in Kinnersley tetrad *)
    dundr  = -( ((2*rp)*En) - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*d\[CapitalSigma]dr + 
             ((rp^2+a^2)*dEndr - a*dLzdr))/(2*\[CapitalSigma]);
    dundx  = -( - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*d\[CapitalSigma]dx + 
             ((rp^2+a^2)*dEndx - a*dLzdx))/(2*\[CapitalSigma]);
    umb1 = ( (-I*(a*(En1-h1[rp,0]) - (Lz1+h2[rp,0]))))/(-Sqrt[2]*\[Zeta]);
    dumbdr = ( - (-I*(a*En - Lz) )*((1-I*a*dzpdr)/\[Zeta]) + 
             (-I*(a*dEndr - dLzdr) + dUzpdr))/(-Sqrt[2]*\[Zeta]);
    dumbdx = ( - (-I*(a*En - Lz))*((-I*a*dzpdx)/\[Zeta]) + 
             (-I*(a*dEndx - dLzdx) + dUzpdx))/(-Sqrt[2]*\[Zeta]);
    Ann0S  = (2*un1)*un + Cdnn; (* Source term *)
    Annt\[Phi]S = (St\[Phi]n)*un;
    AnnrS  = (Srn)*un;
    Ann\[Theta]S  = (S\[Theta]n)*un;
    Anmb0S  = (un1*umb + un*umb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*umb + St\[Phi]mb*un)/2);
    AnmbrS  = ((Srn*umb + Srmb*un)/2);
    Anmb\[Theta]S  = ((S\[Theta]n*umb + S\[Theta]mb*un)/2);
    Ambmb0S  = (2*umb1)*umb + Cdmbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb)*umb;
    AmbmbrS  = (Srmb)*umb;
    Ambmb\[Theta]S  = (S\[Theta]mb)*umb;
    {sumPlus0, sumMinus0, sumPlus1, sumMinus1, dsumPlusdr, dsumMinusdr, dsumPlusdx, dsumMinusdx} = 
      {\[CapitalSigma]*((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RIn - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRIndr + Cmmbmb*fmbmb2*d2RIndr2), (* Total of all quadrants *)
       \[CapitalSigma]*((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUpdr + Cmmbmb*fmbmb2*d2RUpdr2), 
       \[CapitalSigma]*(
        ((Ann0S + Annt\[Phi]S)*RIn*fnn0 + AnnrS*(dRIndr*fnn0 + RIn*dfnn0dr) + Ann\[Theta]S*RIn*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RIn*fnmb0 + AnmbrS*( RIn*dfnmb0dr +  dRIndr*fnmb0) + Anmb\[Theta]S* RIn*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRIndr*fnmb1 + AnmbrS*(dRIndr*dfnmb1dr + d2RIndr2*fnmb1) + Anmb\[Theta]S*dRIndr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RIn*fmbmb0 + AmbmbrS*(  dRIndr*fmbmb0 +   RIn*dfmbmb0dr) + Ambmb\[Theta]S*  RIn*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRIndr*fmbmb1 + AmbmbrS*( d2RIndr2*fmbmb1 +  dRIndr*dfmbmb1dr) + Ambmb\[Theta]S* dRIndr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RIndr2*fmbmb2 + AmbmbrS*(d3RIndr3*fmbmb2 + d2RIndr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RIndr2*dfmbmb2d\[Theta])) + 
        \[CapitalSigma]*((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*dRInd\[Omega] - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*d2RIndrd\[Omega] + Cmmbmb*fmbmb2*d3RIndr2d\[Omega]
          + (Cmnn*dfnn0d\[Omega] + Cmnmb*dfnmb0d\[Omega] + Cmmbmb*dfmbmb0d\[Omega])*RIn - (Cmnmb*dfnmb1d\[Omega] + Cmmbmb*dfmbmb1d\[Omega])*dRIndr + Cmmbmb*dfmbmb2d\[Omega]*d2RIndr2)*\[Omega]1, 
       \[CapitalSigma]*(
        ((Ann0S + Annt\[Phi]S)*RUp*fnn0 + AnnrS*(dRUpdr*fnn0 + RUp*dfnn0dr) + Ann\[Theta]S*RUp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUp*fnmb0 + AnmbrS*( RUp*dfnmb0dr +  dRUpdr*fnmb0) + Anmb\[Theta]S* RUp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUpdr*fnmb1 + AnmbrS*(dRUpdr*dfnmb1dr + d2RUpdr2*fnmb1) + Anmb\[Theta]S*dRUpdr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUp*fmbmb0 + AmbmbrS*(  dRUpdr*fmbmb0 +   RUp*dfmbmb0dr) + Ambmb\[Theta]S*  RUp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUpdr*fmbmb1 + AmbmbrS*( d2RUpdr2*fmbmb1 +  dRUpdr*dfmbmb1dr) + Ambmb\[Theta]S* dRUpdr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RUpdr2*fmbmb2 + AmbmbrS*(d3RUpdr3*fmbmb2 + d2RUpdr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RUpdr2*dfmbmb2d\[Theta])) + 
        \[CapitalSigma]*((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*dRUpd\[Omega] - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*d2RUpdrd\[Omega] + Cmmbmb*fmbmb2*d3RUpdr2d\[Omega]
          + (Cmnn*dfnn0d\[Omega] + Cmnmb*dfnmb0d\[Omega] + Cmmbmb*dfmbmb0d\[Omega])*RUp - (Cmnmb*dfnmb1d\[Omega] + Cmmbmb*dfmbmb1d\[Omega])*dRUpdr + Cmmbmb*dfmbmb2d\[Omega]*d2RUpdr2)*\[Omega]1, 
      \[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Cmnn*fnn0*RIn + Cmnmb*(fnmb0*RIn - fnmb1*dRIndr) + Cmmbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundr*fnn0*RIn + (un*dumbdr+dundr*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdr*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Cmnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RIn + fnn0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)) + 
               Cmnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RIn + fnmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRIndr - fnmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr)) + 
              Cmmbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RIn + fmbmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRIndr - fmbmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RIndr2 + fmbmb2*(d3RIndr3 + d3RIndr2d\[Omega]*d\[Omega]dr)))) 
              ), 
       \[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Cmnn*fnn0*RUp + Cmnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Cmmbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundr*fnn0*RUp + (un*dumbdr+dundr*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdr*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Cmnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RUp + fnn0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)) + 
               Cmnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RUp + fnmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fnmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr)) + 
              Cmmbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RUp + fmbmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fmbmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RUpdr2 + fmbmb2*(d3RUpdr3 + d3RUpdr2d\[Omega]*d\[Omega]dr)))) 
              ), 
      \[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Cmnn*fnn0*RIn + Cmnmb*(fnmb0*RIn - fnmb1*dRIndr) + Cmmbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundx*fnn0*RIn + (un*dumbdx+dundx*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdx*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Cmnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RIn + fnn0*(dRInd\[Omega]*d\[Omega]dx)) + 
               Cmnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RIn + fnmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRIndr - fnmb1*(d2RIndrd\[Omega]*d\[Omega]dx)) + 
              Cmmbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RIn + fmbmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfmbmb1d\[Theta]*d\[Theta]pdx + dfmbmb1d\[Omega]*d\[Omega]dx)*dRIndr - fmbmb1*(d2RIndrd\[Omega]*d\[Omega]dx) + 
                      (dfmbmb2d\[Theta]*d\[Theta]pdx + dfmbmb2d\[Omega]*d\[Omega]dx)*d2RIndr2 + fmbmb2*(d3RIndr2d\[Omega]*d\[Omega]dx)))) 
              ), 
       \[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Cmnn*fnn0*RUp + Cmnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Cmmbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundx*fnn0*RUp + (un*dumbdx+dundx*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdx*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Cmnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RUp + fnn0*(dRUpd\[Omega]*d\[Omega]dx)) + 
               Cmnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RUp + fnmb0*(dRUpd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRUpdr - fnmb1*(d2RUpdrd\[Omega]*d\[Omega]dx)) + 
              Cmmbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RUp + fmbmb0*(dRUpd\[Omega]*d\[Omega]dx)
                    - (dfmbmb1d\[Theta]*d\[Theta]pdx + dfmbmb1d\[Omega]*d\[Omega]dx)*dRUpdr - fmbmb1*(d2RUpdrd\[Omega]*d\[Omega]dx) + 
                      (dfmbmb2d\[Theta]*d\[Theta]pdx + dfmbmb2d\[Omega]*d\[Omega]dx)*d2RUpdr2 + fmbmb2*(d3RUpdr2d\[Omega]*d\[Omega]dx)))) 
              )
       };
  W = (RIn*dRUpdr - dRIndr*RUp)/\[CapitalDelta]; (* Invariant Wronskian *)
  dWd\[Omega]= (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/\[CapitalDelta]; (* \[Omega]-derivative of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*\[Omega]1/W*sumPlus0 )/(\[CapitalGamma]*W); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*\[Omega]1/W*sumMinus0)/(\[CapitalGamma]*W);
  dCPlusdr  = 2*Pi*(dsumPlusdr  - d\[CapitalGamma]dr/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*d\[Omega]dr/W*sumPlus0 )/(\[CapitalGamma]*W); (* Derivatives of the amplitudes *)
  dCMinusdr = 2*Pi*(dsumMinusdr - d\[CapitalGamma]dr/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*d\[Omega]dr/W*sumMinus0)/(\[CapitalGamma]*W);
  dCPlusdx  = 2*Pi*(dsumPlusdx  - d\[CapitalGamma]dx/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*d\[Omega]dx/W*sumPlus0 )/(\[CapitalGamma]*W);
  dCMinusdx = 2*Pi*(dsumMinusdx - d\[CapitalGamma]dx/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*d\[Omega]dx/W*sumMinus0)/(\[CapitalGamma]*W);
  <|
    "l" -> l,
    "m" -> m,
    "k" -> 0,
    "n" -> 0,
    "\[Omega]" -> \[Omega],
    "Amplitudes" -> <|
      "\[ScriptCapitalI]" -> CPlus0,
      "\[ScriptCapitalH]" -> CMinus0
    |>,
    "AmplitudesCorrection" -> <|
      "\[ScriptCapitalI]" -> CPlus1,
      "\[ScriptCapitalH]" -> CMinus1
    |>,
    "AmplitudesDerivatives" -> <|
      "r" -> <|
        "\[ScriptCapitalI]" -> dCPlusdr,
        "\[ScriptCapitalH]" -> dCMinusdr
      |>,
      "x" -> <|
        "\[ScriptCapitalI]" -> dCPlusdx,
        "\[ScriptCapitalH]" -> dCMinusdx
      |>
    |>,
    "\[Alpha]" -> \[Alpha],
    "S" -> N[SWSH[Pi/2,0][[1]]],
    "Fluxes" -> <|
      "Energy" -> <|
        "\[ScriptCapitalI]" -> Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]" -> \[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum" -> <|
        "\[ScriptCapitalI]" -> Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]" -> \[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection" -> <|
      "Energy" -> <|
        "\[ScriptCapitalI]" -> (2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]" -> \[Alpha]*(2*Re[CMinus1*Conjugate[CMinus0]] + (d\[Alpha]d\[Omega]/\[Alpha]*Abs[CMinus0]^2 - 2*Abs[CMinus0]^2/\[Omega])*\[Omega]1)/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum" -> <|
        "\[ScriptCapitalI]" -> (2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesDerivatives" -> <|
      "r" -> <|
        "Energy"  ->  <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdr*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*d\[Omega]dr/\[Omega])/(4Pi*\[Omega]^2),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dr/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdr*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*d\[Omega]dr/\[Omega])/(4Pi*\[Omega]^2)
        |>,
        "AngularMomentum" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdr*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*d\[Omega]dr/\[Omega])*m/(4Pi*\[Omega]^3),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dr/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdr*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*d\[Omega]dr/\[Omega])*m/(4Pi*\[Omega]^3)
        |>
      |>,
      "x" -> <|
        "Energy" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdx*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*d\[Omega]dx/\[Omega])/(4Pi*\[Omega]^2),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dx/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdx*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*d\[Omega]dx/\[Omega])/(4Pi*\[Omega]^2)
        |>,
        "AngularMomentum" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdx*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*d\[Omega]dx/\[Omega])*m/(4Pi*\[Omega]^3),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dx/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdx*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*d\[Omega]dx/\[Omega])*m/(4Pi*\[Omega]^3)
        |>
      |>
    |>
  |>
]*)


(* ::Subsection::Closed:: *)
(*Equatorial*)


Options[TeukolskySpinModeEquatorialCorrectionAnalytical] = {WorkingPrecision->30};
TeukolskySpinModeEquatorialCorrectionAnalytical[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,orbit_,{angparNew_,RCorrection_},OptionsPattern[]]:=Module[{a,p,e,x,
    En0,Lz0,K0,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],r,\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t1,\[Omega],\[Omega]1,SWSH,dSWSHd\[Omega],R,\[Lambda],d\[Lambda]d\[Omega],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumPlus1,sumMinus0,sumMinus1,stepsr,ir,wr,
    rp,Urp,expr,expr1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],V,dVd\[Omega],
    Rp,RIn,dRIndr,d2RIndr2,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,dDRInd\[Omega],d2DRIndrd\[Omega],dDDRInd\[Omega],dDRUpd\[Omega],d2DRUpdrd\[Omega],dDDRUpd\[Omega],
    zp,\[Theta]p,sin\[Theta]p,K\[Theta],dK\[Theta]d\[Omega],S,dSd\[Theta],d2Sd\[Theta]2,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],
    L2S,dL2Sd\[Theta],L1L2S,dL2Sd\[Omega],d2L2Sd\[Theta]d\[Omega],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,vn,vmb,
    FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega],
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega],
    CPlus0,CPlus1,CMinus0,CMinus1,
    \[ScriptCapitalF]En\[ScriptCapitalI],\[ScriptCapitalF]En\[ScriptCapitalH],\[ScriptCapitalF]Lz\[ScriptCapitalI],\[ScriptCapitalF]Lz\[ScriptCapitalH],\[ScriptCapitalF]En\[ScriptCapitalI]1,\[ScriptCapitalF]En\[ScriptCapitalH]1,\[ScriptCapitalF]Lz\[ScriptCapitalI]1,\[ScriptCapitalF]Lz\[ScriptCapitalH]1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["Inclination"];
  If[x!=1,Return[$Failed]];
  En0 = orbit["Energy"]; (* Shifted constants of motion *)
  Lz0 = orbit["AngularMomentum"];
  K0 = orbit["CarterConstant"] + (Lz0-a*En0)^2;
  {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[orbit["Frequencies"]];
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t; (* BL frequencies *)
  r = orbit["Trajectory"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  \[CapitalUpsilon]\[Tau] = orbit["ProperTimeFrequency"];
  \[CapitalUpsilon]t1 = -3*\[CapitalUpsilon]\[Tau]/(2*Sqrt[K0]);
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r; (* Frequency of mode *)
  \[Omega]1 = 3*\[CapitalUpsilon]\[Tau]*\[Omega]/(2*Sqrt[K0]*\[CapitalUpsilon]t);
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];  
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,
                                     SetPrecision[a, OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega], OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                         "precODE" -> OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  Print["Calculated R: "<>ToString[AbsoluteTiming[
    R = RCorrection[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],
                         SetPrecision[\[Omega],OptionValue[WorkingPrecision]+5],1,
                         SetPrecision[\[Lambda],OptionValue[WorkingPrecision]+5],
                         SetPrecision[d\[Lambda]d\[Omega],OptionValue[WorkingPrecision]+5],e,p,"precODE"->OptionValue[WorkingPrecision]]
                         ][[1]]]];
  \[Lambda]1 = d\[Lambda]d\[Omega]*\[Omega]1;
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  zp = 0;
  \[Theta]p = N[Pi/2,Precision[{a,p,e}]];
  sin\[Theta]p = 1;
  K\[Theta] = a*\[Omega]*(1-zp^2) - m;
  dK\[Theta]d\[Omega] = a*(1-zp^2);
  (* Spin-weighted spheroidal harmonics and their \[Theta]-derivatives *)
  {S,dSd\[Theta]}=SWSH[\[Theta]p]; 
  d2Sd\[Theta]2 = -(zp/sin\[Theta]p)*dSd\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
  (* \[Omega]-derivatives of S *)
  {dSd\[Omega],d2Sd\[Theta]d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
  d3Sd\[Theta]2d\[Omega] = -(zp/sin\[Theta]p)*d2Sd\[Theta]d\[Omega] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega]*(1-zp^2) + 4*a*zp + 2*m*a + d\[Lambda]d\[Omega])*S;
  (* Operators L and their \[Omega]-derivatives *)
  L2S = dSd\[Theta] + (K\[Theta] + 2 zp)*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
  dL2Sd\[Theta] = d2Sd\[Theta]2 + (K\[Theta] + 2 zp)*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*S;
  L1L2S = dL2Sd\[Theta] + (K\[Theta] + zp)*L2S/sin\[Theta]p; 
  dL2Sd\[Omega] = d2Sd\[Theta]d\[Omega] + (K\[Theta] + 2 zp)*dSd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
  d2L2Sd\[Theta]d\[Omega] = d3Sd\[Theta]2d\[Omega] + (K\[Theta] + 2 zp)*d2Sd\[Theta]d\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*dSd\[Omega] + (2*a*zp - dK\[Theta]d\[Omega]/sin\[Theta]p^2*zp)*S;
  dL1L2Sd\[Omega] = d2L2Sd\[Theta]d\[Omega] + (K\[Theta] + zp)*dL2Sd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*L2S/sin\[Theta]p;
  (* numbers of steps for wr and wz integration *)
  stepsr = Max[2^(5+Ceiling[Log2[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]]]),
               2^(5+Ceiling[Log2[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]]]),32];
  Print[ToString[stepsr]<>" steps in wr"];
  {sumPlus0, sumPlus1, sumMinus0, sumMinus1} = Sum[(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e}]];
    rp = r[wr];
    Urp = {1,-1}*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+K0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    expr1 = {1,-1}*I*(\[Omega]1*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr] + \[Omega]*(-3/(2*Sqrt[K0]))*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]r"][wr]);(* Linear part of the exponential term *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    dKdr = 2*rp*\[Omega];
    dKd\[Omega]  = (rp^2+a^2);
    d2Kdrd\[Omega] = 2*rp;
    d\[CapitalDelta] = 2*(rp-1);
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
    (* In solutions of the radial equation and their r and \[Omega] derivatives *)
    Rp = R[rp];
    RIn    = Rp["R0"]["In"]["R"]; 
    dRIndr   = Rp["R0"]["In"]["dR"];
    d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dRInd\[Omega]   = Rp["R1"]["In"]["R"]; (* \[Omega]-derivatives of the radial function *)
    d2RIndrd\[Omega]  = Rp["R1"]["In"]["dR"];
    d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
    (* Operators D and their r and \[Omega] derivatives *)
    DRIn = dRIndr - I*K/\[CapitalDelta]*RIn;
    dDRInd\[Omega] = d2RIndrd\[Omega] - I*K/\[CapitalDelta]*dRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RIn;
    dDRIndr = d2RIndr2 - I*K/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
    d2DRIndrd\[Omega] = d3RIndr2d\[Omega] - I*K/\[CapitalDelta]*d2RIndrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRInd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
    DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
    dDDRInd\[Omega] = d2DRIndrd\[Omega] - I*K/\[CapitalDelta]*dDRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRIn;
    (* Up solutions of the radial equation and their r and \[Omega] derivatives *)
    RUp    = Rp["R0"]["Up"]["R"]; 
    dRUpdr   = Rp["R0"]["Up"]["dR"];
    d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dRUpd\[Omega]   = Rp["R1"]["Up"]["R"]; (* \[Omega]-derivative of the radial function *)
    d2RUpdrd\[Omega]  = Rp["R1"]["Up"]["dR"];
    d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
    (* Operators D and their r and \[Omega] derivatives *)
    DRUp = dRUpdr - I*K/\[CapitalDelta]*RUp;
    dDRUpd\[Omega] = d2RUpdrd\[Omega] - I*K/\[CapitalDelta]*dRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RUp;
    dDRUpdr = d2RUpdr2 - I*K/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
    d2DRUpdrd\[Omega] = d3RUpdr2d\[Omega] - I*K/\[CapitalDelta]*d2RUpdrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRUpd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
    DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
    dDDRUpd\[Omega] = d2DRUpdrd\[Omega] - I*K/\[CapitalDelta]*dDRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRUp;
    (* Integration f functions of r and \[Theta] *)
    \[Zeta] = rp-I*a*zp;
    \[Zeta]bar = rp+I*a*zp;
    vn = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
    vmb = (-I*(a*sin\[Theta]p^2*En0 - Lz0))/sin\[Theta]p^2;
    (* In solution *)
    {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = FabEq[rp,a,RIn,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
    {dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RIn,DRIn,DDRIn,dRInd\[Omega],dDRInd\[Omega],dDDRInd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
    (* Up solution *)
    {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = FabEq[rp,a,RUp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
    {dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RUp,DRUp,DDRUp,dRUpd\[Omega],dDRUpd\[Omega],dDDRUpd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
    {Total[(vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr^{1,-1}],
    Total[(
        (vn*GnIn + vmb*GmbIn)/Sqrt[K0] +
        (vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr1 +
        (vn*vn*dFnnInd\[Omega] + vn*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
     )*expr^{1,-1}],
     Total[(vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr^{1,-1}],
     Total[(
        (vn*GnUp + vmb*GmbUp)/Sqrt[K0] +
        (vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr1 +
        (vn*vn*dFnnUpd\[Omega] + vn*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
       )*expr^{1,-1}]},
       {ir, 1, stepsr/2}
  ];
  W = (RIn*dRUpdr - dRIndr*RUp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  W1 = (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/(rp^2-2*rp+a^2)*\[Omega]1; (* Invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalUpsilon]t*W*stepsr); (* Amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalUpsilon]t*W*stepsr);
  CPlus1  = 2*Pi*sumPlus1/(\[CapitalUpsilon]t*W*stepsr) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CPlus0; (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*sumMinus1/(\[CapitalUpsilon]t*W*stepsr) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CMinus0;
  \[ScriptCapitalF]En\[ScriptCapitalI] = Abs[CPlus0]^2/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]En\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI] = Abs[CPlus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]En\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]En\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  <|
    "l"->l, "m"->m, "k"->0, "n"->n,
    "\[Omega]"->\[Omega], "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"-><|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"-><|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha], "\[Alpha]Correction"->\[Alpha]1,
    "S"->S,
    "Fluxes"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]|>
    |>,
    "FluxesCorrection"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]1|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]1|>
    |>,
    "stepsr"->stepsr
  |>
]


Trajectory[orbit_]:=Module[{M=1,a,En,L,Q,r1,r2,r3,r4,kr,rp,rm,hr,hp,hm,traj},
  a = orbit["a"];
  {En,L,Q} = Values[orbit["ConstantsOfMotion"]];
  {r1,r2,r3,r4} = orbit["RadialRoots"];	
  kr = (r1-r2)*(r3-r4)/((r1-r3)*(r2-r4));
  rp=M+Sqrt[M^2-a^2];
  rm=M-Sqrt[M^2-a^2];
  hr=(r1-r2)/(r1-r3);
  hp=((r1-r2)(r3-rp))/((r1-r3)(r2-rp));
  hm=((r1-r2)(r3-rm))/((r1-r3)(r2-rm));
  
  traj = Function[{qr}, Module[{sn,\[Psi]r,\[CapitalDelta]E,\[CapitalDelta]\[CapitalPi]hr,\[CapitalDelta]\[CapitalPi]hm,\[CapitalDelta]\[CapitalPi]hp,r,Ur,\[CapitalDelta]tr,\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Tau]r},
    sn = JacobiSN[EllipticK[kr]/\[Pi] qr,kr];
    (*cn = JacobiCN[EllipticK[kr]/\[Pi] qr,kr];
    dn = JacobiDN[EllipticK[kr]/\[Pi] qr,kr];*)
    \[Psi]r = ArcSin[sn];
    \[CapitalDelta]E   = EllipticE[kr] qr/\[Pi]-EllipticE[\[Psi]r,kr];
    \[CapitalDelta]\[CapitalPi]hr = EllipticPi[hr,kr] qr/\[Pi]-EllipticPi[hr,\[Psi]r,kr];
    \[CapitalDelta]\[CapitalPi]hm = EllipticPi[hm,kr] qr/\[Pi]-EllipticPi[hm,\[Psi]r,kr];
    \[CapitalDelta]\[CapitalPi]hp = EllipticPi[hp,kr] qr/\[Pi]-EllipticPi[hp,\[Psi]r,kr];
    
    r = (r3(r1 - r2)sn^2-r2(r1-r3))/((r1-r2)sn^2-(r1-r3));
    Ur = Sqrt[(1-En^2)/(r1-r3)*(r2-r4)]*(r1-r2) (r2-r3) Sqrt[(1-sn^2)*(1-kr*sn^2)] sn/(1-hr sn^2)^2;
    \[CapitalDelta]tr = -1/Sqrt[(1-En^2) (r1-r3) (r2-r4)] (
      - 4 (r2-r3)/(rp-rm) (
      -(-2 a^2 En + rm (4 En - a L))/((-rm+r2) (-rm+r3)) \[CapitalDelta]\[CapitalPi]hm 
      +(-2 a^2 En + rp (4 En - a L))/((-rp+r2) (-rp+r3)) \[CapitalDelta]\[CapitalPi]hp) 
      + 2 En (r2-r3) (2 + 1/(1-En^2)) \[CapitalDelta]\[CapitalPi]hr 
      + En (r2-r4) ((r1-r3) \[CapitalDelta]E + (r1-r2) sn/(1-hr sn^2) Sqrt[(1-sn^2) (1-kr sn^2)]) );

    \[CapitalDelta]\[Phi]r = (2 a (r2-r3) (
      -(2 rm En - a L)/((-rm+r2) (-rm+r3)) \[CapitalDelta]\[CapitalPi]hm
      +(2 rp En - a L)/((-rp+r2) (-rp+r3)) \[CapitalDelta]\[CapitalPi]hp)
      )/((-rm+rp) Sqrt[(1-En^2) (r1-r3) (r2-r4)]);

    \[CapitalDelta]\[Tau]r = -1/Sqrt[(1-En^2) (r1-r3) (r2-r4)] ((r2-r3) 2/(1-En^2) \[CapitalDelta]\[CapitalPi]hr
      +(r2-r4) ((r1-r3) \[CapitalDelta]E + (r1-r2) sn/(1-hr sn^2) Sqrt[(1-sn^2) (1-kr sn^2)]) );
    <|"r"->r, "Ur"->Ur, "\[CapitalDelta]tr"->\[CapitalDelta]tr, "\[CapitalDelta]\[Phi]r"->\[CapitalDelta]\[Phi]r, "\[CapitalDelta]\[Tau]r"->\[CapitalDelta]\[Tau]r|>
  ]];
  <|
    "a"->a,
    "p"->orbit["p"],
    "e"->orbit["e"],
    "Inclination"->orbit["Inclination"],
    "Energy"->En,
    "AngularMomentum"->L,
    "CarterConstant"->Q,
    "Frequencies"->orbit["Frequencies"],
    "ProperTimeFrequency"->orbit["ProperTimeFrequency"],
    "Trajectory"->traj,
    "TrajectoryDeltas"->orbit["TrajectoryDeltas"]
  |>
]


Options[TeukolskySpinModeEquatorialCorrectionAnalyticalNew] = {WorkingPrecision->30};
TeukolskySpinModeEquatorialCorrectionAnalyticalNew[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,orbit_,{angparNew_,RCorrection_},OptionsPattern[]]:=Module[{prec,a,p,e,x,
    En0,Lz0,K0,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],traj,\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t1,\[Omega],\[Omega]1,SWSH,dSWSHd\[Omega],R,Rp,\[Lambda],d\[Lambda]d\[Omega],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumPlus1,sumMinus0,sumMinus1,stepsr,wr,
    trajp,rp,Urp,expr,expr1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],V,dVd\[Omega],
    RIn,dRIndr,d2RIndr2,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,dDRInd\[Omega],d2DRIndrd\[Omega],dDDRInd\[Omega],dDRUpd\[Omega],d2DRUpdrd\[Omega],dDDRUpd\[Omega],
    \[Theta]p,K\[Theta],dK\[Theta]d\[Omega],S,dSd\[Theta],d2Sd\[Theta]2,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],
    L2S,dL2Sd\[Theta],L1L2S,dL2Sd\[Omega],d2L2Sd\[Theta]d\[Omega],dL1L2Sd\[Omega],
    vnplus,vnminus,vmb,
    FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega],
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega],
    CPlus0,CPlus1,CMinus0,CMinus1,
    \[ScriptCapitalF]En\[ScriptCapitalI],\[ScriptCapitalF]En\[ScriptCapitalH],\[ScriptCapitalF]Lz\[ScriptCapitalI],\[ScriptCapitalF]Lz\[ScriptCapitalH],\[ScriptCapitalF]En\[ScriptCapitalI]1,\[ScriptCapitalF]En\[ScriptCapitalH]1,\[ScriptCapitalF]Lz\[ScriptCapitalI]1,\[ScriptCapitalF]Lz\[ScriptCapitalH]1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  prec = OptionValue[WorkingPrecision];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["Inclination"];
  If[x!=1,Return[$Failed]];
  En0 = orbit["Energy"]; (* Shifted constants of motion *)
  Lz0 = orbit["AngularMomentum"];
  K0 = orbit["CarterConstant"] + (Lz0-a*En0)^2;
  {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[orbit["Frequencies"]];(* Mino frequencies *)
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t; (* BL frequencies *)
  traj = orbit["Trajectory"];(* Function for the trajectory *)
  \[CapitalUpsilon]\[Tau] = orbit["ProperTimeFrequency"];
  \[CapitalUpsilon]t1 = -3*\[CapitalUpsilon]\[Tau]/(2*Sqrt[K0]);
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r; (* Frequency of mode *)
  \[Omega]1 = 3*\[CapitalUpsilon]\[Tau]*\[Omega]/(2*Sqrt[K0]*\[CapitalUpsilon]t);(* Linear part of the frequency *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,
                                     SetPrecision[a, prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega], prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                         "precODE" -> prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  Print["Calculated R: "<>ToString[AbsoluteTiming[
    R = RCorrection[-2,l,m,SetPrecision[a,prec+5],
                         SetPrecision[\[Omega],prec+5],1,
                         SetPrecision[\[Lambda],prec+5],
                         SetPrecision[d\[Lambda]d\[Omega],prec+5],e,p,"precODE"->prec];
  ][[1]]]];
  \[Lambda]1 = d\[Lambda]d\[Omega]*\[Omega]1;
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  (* Functions of \[Theta] for equatorial orbits *)
  \[Theta]p = N[Pi/2,prec];
  K\[Theta] = a*\[Omega] - m;
  dK\[Theta]d\[Omega] = a;
  (* Spin-weighted spheroidal harmonics and their \[Theta]-derivatives *)
  {S,dSd\[Theta]}=SWSH[\[Theta]p]; 
  d2Sd\[Theta]2 = - (-a^2*\[Omega]^2 - m^2 - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
  (* \[Omega]-derivatives of S *)
  {dSd\[Omega],d2Sd\[Theta]d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta])  *)
  d3Sd\[Theta]2d\[Omega] = - (-a^2*\[Omega]^2 - m^2 - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega] + 2*m*a + d\[Lambda]d\[Omega])*S;
  (* Operators L and their \[Omega]-derivatives *)
  L2S = dSd\[Theta] + K\[Theta]*S;
  dL2Sd\[Theta] = d2Sd\[Theta]2 + K\[Theta]*dSd\[Theta] -2*S;
  L1L2S = dL2Sd\[Theta] + K\[Theta]*L2S; 
  dL2Sd\[Omega] = d2Sd\[Theta]d\[Omega] + K\[Theta]*dSd\[Omega] + dK\[Theta]d\[Omega]*S;
  d2L2Sd\[Theta]d\[Omega] = d3Sd\[Theta]2d\[Omega] + K\[Theta]*d2Sd\[Theta]d\[Omega] + dK\[Theta]d\[Omega]*dSd\[Theta] -2*dSd\[Omega];
  dL1L2Sd\[Omega] = d2L2Sd\[Theta]d\[Omega] + K\[Theta]*dL2Sd\[Omega] + dK\[Theta]d\[Omega]*L2S;
  (* numbers of steps for wr integration from the exponential functions *)
  stepsr = Max[2^(5+Ceiling[Log2[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi]+n)]]]),
               2^(5+Ceiling[Log2[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0 ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0 ]+n)]]]),32];
  Print[ToString[stepsr]<>" steps in wr"];
  wr = N[Range[Pi/stepsr,(stepsr-1)*Pi/stepsr,2Pi/stepsr],OptionValue[WorkingPrecision]];
  Print["traj: "<>ToString[AbsoluteTiming[
    trajp = traj[wr];
  ][[1]]]];
  rp = trajp["r"];
  Urp = trajp["Ur"];(* Geodesic radial velocity *)
  expr = Exp[I*(\[Omega]*(trajp["\[CapitalDelta]tr"])-m*(trajp["\[CapitalDelta]\[Phi]r"])+n*wr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
  expr1 = I*(\[Omega]1*trajp["\[CapitalDelta]tr"] + \[Omega]*(-3/(2*Sqrt[K0]))*trajp["\[CapitalDelta]\[Tau]r"]);(* Linear part of the exponential term *)
  Print["Rp: "<>ToString[AbsoluteTiming[
    Rp = R[rp];(* Evaluate the radial functions and their derivatives *)
  ][[1]]]];
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKdr = 2*rp*\[Omega];
  dKd\[Omega]  = (rp^2+a^2);
  d2Kdrd\[Omega] = 2*rp;
  d\[CapitalDelta] = 2*(rp-1);
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  (* In solutions of the radial equation and their r and \[Omega] derivatives *)
  RIn    = Rp["R0"]["In"]["R"]; 
  dRIndr   = Rp["R0"]["In"]["dR"];
  d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  dRInd\[Omega]   = Rp["R1"]["In"]["R"]; (* \[Omega]-derivatives of the radial function *)
  d2RIndrd\[Omega]  = Rp["R1"]["In"]["dR"];
  d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRIn = dRIndr - I*K/\[CapitalDelta]*RIn;
  dDRInd\[Omega] = d2RIndrd\[Omega] - I*K/\[CapitalDelta]*dRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RIn;
  dDRIndr = d2RIndr2 - I*K/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  d2DRIndrd\[Omega] = d3RIndr2d\[Omega] - I*K/\[CapitalDelta]*d2RIndrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRInd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
  dDDRInd\[Omega] = d2DRIndrd\[Omega] - I*K/\[CapitalDelta]*dDRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRIn;
  (* Up solutions of the radial equation and their r and \[Omega] derivatives *)
  RUp    = Rp["R0"]["Up"]["R"]; 
  dRUpdr   = Rp["R0"]["Up"]["dR"];
  d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  dRUpd\[Omega]   = Rp["R1"]["Up"]["R"]; (* \[Omega]-derivative of the radial function *)
  d2RUpdrd\[Omega]  = Rp["R1"]["Up"]["dR"];
  d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRUp = dRUpdr - I*K/\[CapitalDelta]*RUp;
  dDRUpd\[Omega] = d2RUpdrd\[Omega] - I*K/\[CapitalDelta]*dRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RUp;
  dDRUpdr = d2RUpdr2 - I*K/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  d2DRUpdrd\[Omega] = d3RUpdr2d\[Omega] - I*K/\[CapitalDelta]*d2RUpdrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRUpd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
  dDDRUpd\[Omega] = d2DRUpdrd\[Omega] - I*K/\[CapitalDelta]*dDRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRUp;
  (* components of the velocity *)
  vnplus = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
  vnminus = -((rp^2+a^2)*En0 - a*Lz0 - Urp)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
  vmb = (-I*(a*En0 - Lz0));
  (* In solution *)
  {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = FabEq[rp,a,RIn,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
  {dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RIn,DRIn,DDRIn,dRInd\[Omega],dDRInd\[Omega],dDDRInd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
  (* Up solution *)
  {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = FabEq[rp,a,RUp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
  {dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RUp,DRUp,DDRUp,dRUpd\[Omega],dDRUpd\[Omega],dDDRUpd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
  sumPlus0 = Total[(vnplus*vnplus*FnnIn + vnplus*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr + (vnminus*vnminus*FnnIn + vnminus*vmb*FnmbIn + vmb*vmb*FmbmbIn)/expr,"Method"->"CompensatedSummation"];
  sumPlus1 = Total[(
      (vnplus*GnIn + vmb*GmbIn)/Sqrt[K0] +
      (vnplus*vnplus*FnnIn + vnplus*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr1 +
      (vnplus*vnplus*dFnnInd\[Omega] + vnplus*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
   )*expr + 
   (
      (vnminus*GnIn + vmb*GmbIn)/Sqrt[K0] +
      (vnminus*vnminus*FnnIn + vnminus*vmb*FnmbIn + vmb*vmb*FmbmbIn)*(-expr1) +
      (vnminus*vnminus*dFnnInd\[Omega] + vnminus*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
   )/expr,"Method"->"CompensatedSummation"];
  sumMinus0 = Total[(vnplus*vnplus*FnnUp + vnplus*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr + (vnminus*vnminus*FnnUp + vnminus*vmb*FnmbUp + vmb*vmb*FmbmbUp)/expr,"Method"->"CompensatedSummation"];
  sumMinus1 = Total[(
      (vnplus*GnUp + vmb*GmbUp)/Sqrt[K0] +
      (vnplus*vnplus*FnnUp + vnplus*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr1 +
      (vnplus*vnplus*dFnnUpd\[Omega] + vnplus*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
   )*expr + 
   (
      (vnminus*GnUp + vmb*GmbUp)/Sqrt[K0] +
      (vnminus*vnminus*FnnUp + vnminus*vmb*FnmbUp + vmb*vmb*FmbmbUp)*(-expr1) +
      (vnminus*vnminus*dFnnUpd\[Omega] + vnminus*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
   )/expr,"Method"->"CompensatedSummation"];
  W = (RIn[[1]]*dRUpdr[[1]] - dRIndr[[1]]*RUp[[1]])/(rp[[1]]^2-2*rp[[1]]+a^2); (* Invariant Wronskian *)
  W1 = (dRInd\[Omega][[1]]*dRUpdr[[1]] + RIn[[1]]*d2RUpdrd\[Omega][[1]] - d2RIndrd\[Omega][[1]]*RUp[[1]] - dRIndr[[1]]*dRUpd\[Omega][[1]])/(rp[[1]]^2-2*rp[[1]]+a^2)*\[Omega]1; (* derivative of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalUpsilon]t*W*stepsr); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalUpsilon]t*W*stepsr);
  CPlus1  = 2*Pi*sumPlus1/(\[CapitalUpsilon]t*W*stepsr) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CPlus0; (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*sumMinus1/(\[CapitalUpsilon]t*W*stepsr) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CMinus0;
  \[ScriptCapitalF]En\[ScriptCapitalI] = Abs[CPlus0]^2/(4Pi*\[Omega]^2); (* Fluxes and their linear parts *)
  \[ScriptCapitalF]En\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI] = Abs[CPlus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]En\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]En\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  <|
    "l"->l, "m"->m, "k"->0, "n"->n,
    "\[Omega]"->\[Omega], "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"-><|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"-><|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha], "\[Alpha]Correction"->\[Alpha]1,
    "S"->S,
    "Fluxes"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]|>
    |>,
    "FluxesCorrection"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]1|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]1|>
    |>,
    "stepsr"->stepsr
  |>
]


Options[TeukolskySpinModeEquatorialCorrectionAnalyticalNew2] = {WorkingPrecision->30};
TeukolskySpinModeEquatorialCorrectionAnalyticalNew2[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,orbit_,{angparNew_,RCorrection_},OptionsPattern[]]:=Module[{prec,a,p,e,x,
    En0,Lz0,K0,r3,r4,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t1,\[Omega],\[Omega]1,SWSH,dSWSHd\[Omega],R,Rp,\[Lambda],d\[Lambda]d\[Omega],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumPlus1,sumMinus0,sumMinus1,stepsr,\[Chi]r,
    cos\[Chi]r,sin\[Chi]r,rp,d\[Lambda]d\[Chi]r,Urp,\[CapitalDelta]tr,\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Tau]r,expr,expr1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],V,dVd\[Omega],
    RIn,dRIndr,d2RIndr2,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,dDRInd\[Omega],d2DRIndrd\[Omega],dDDRInd\[Omega],dDRUpd\[Omega],d2DRUpdrd\[Omega],dDDRUpd\[Omega],
    \[Theta]p,K\[Theta],dK\[Theta]d\[Omega],S,dSd\[Theta],d2Sd\[Theta]2,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],
    L2S,dL2Sd\[Theta],L1L2S,dL2Sd\[Omega],d2L2Sd\[Theta]d\[Omega],dL1L2Sd\[Omega],
    vnplus,vnminus,vmb,
    FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega],
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega],
    CPlus0,CPlus1,CMinus0,CMinus1,
    \[ScriptCapitalF]En\[ScriptCapitalI],\[ScriptCapitalF]En\[ScriptCapitalH],\[ScriptCapitalF]Lz\[ScriptCapitalI],\[ScriptCapitalF]Lz\[ScriptCapitalH],\[ScriptCapitalF]En\[ScriptCapitalI]1,\[ScriptCapitalF]En\[ScriptCapitalH]1,\[ScriptCapitalF]Lz\[ScriptCapitalI]1,\[ScriptCapitalF]Lz\[ScriptCapitalH]1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  prec = OptionValue[WorkingPrecision];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["Inclination"];
  If[x!=1,Return[$Failed]];
  En0 = orbit["Energy"]; (* Shifted constants of motion *)
  Lz0 = orbit["AngularMomentum"];
  K0 = orbit["CarterConstant"] + (Lz0-a*En0)^2;
  {r3,r4} = orbit["RadialRoots"][[3;;4]];
  {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[orbit["Frequencies"]];(* Mino frequencies *)
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t; (* BL frequencies *)
  \[CapitalUpsilon]\[Tau] = orbit["ProperTimeFrequency"];
  \[CapitalUpsilon]t1 = -3*\[CapitalUpsilon]\[Tau]/(2*Sqrt[K0]);
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r; (* Frequency of mode *)
  \[Omega]1 = 3*\[CapitalUpsilon]\[Tau]*\[Omega]/(2*Sqrt[K0]*\[CapitalUpsilon]t);(* Linear part of the frequency *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,
                                     SetPrecision[a, prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega], prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                         "precODE" -> prec+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  Print["Calculated R: "<>ToString[AbsoluteTiming[
    R = RCorrection[-2,l,m,SetPrecision[a,prec+5],
                         SetPrecision[\[Omega],prec+5],1,
                         SetPrecision[\[Lambda],prec+5],
                         SetPrecision[d\[Lambda]d\[Omega],prec+5],e,p,"precODE"->prec];
  ][[1]]]];
  \[Lambda]1 = d\[Lambda]d\[Omega]*\[Omega]1;
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  (* Functions of \[Theta] for equatorial orbits *)
  \[Theta]p = N[Pi/2,prec];
  K\[Theta] = a*\[Omega] - m;
  dK\[Theta]d\[Omega] = a;
  (* Spin-weighted spheroidal harmonics and their \[Theta]-derivatives *)
  {S,dSd\[Theta]}=SWSH[\[Theta]p]; 
  d2Sd\[Theta]2 = - (-a^2*\[Omega]^2 - m^2 - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
  (* \[Omega]-derivatives of S *)
  {dSd\[Omega],d2Sd\[Theta]d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta])  *)
  d3Sd\[Theta]2d\[Omega] = - (-a^2*\[Omega]^2 - m^2 - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega] + 2*m*a + d\[Lambda]d\[Omega])*S;
  (* Operators L and their \[Omega]-derivatives *)
  L2S = dSd\[Theta] + K\[Theta]*S;
  dL2Sd\[Theta] = d2Sd\[Theta]2 + K\[Theta]*dSd\[Theta] -2*S;
  L1L2S = dL2Sd\[Theta] + K\[Theta]*L2S; 
  dL2Sd\[Omega] = d2Sd\[Theta]d\[Omega] + K\[Theta]*dSd\[Omega] + dK\[Theta]d\[Omega]*S;
  d2L2Sd\[Theta]d\[Omega] = d3Sd\[Theta]2d\[Omega] + K\[Theta]*d2Sd\[Theta]d\[Omega] + dK\[Theta]d\[Omega]*dSd\[Theta] -2*dSd\[Omega];
  dL1L2Sd\[Omega] = d2L2Sd\[Theta]d\[Omega] + K\[Theta]*dL2Sd\[Omega] + dK\[Theta]d\[Omega]*L2S;
  (* numbers of steps for wr integration from the exponential functions *)
  stepsr = Max[2^(5+Ceiling[Log2[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi]+n)]]]),
               2^(5+Ceiling[Log2[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0 ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0 ]+n)]]]),32];
  Print[ToString[stepsr]<>" steps in wr"];
  \[Chi]r = N[Range[Pi/stepsr,(stepsr-1)*Pi/stepsr,2Pi/stepsr],OptionValue[WorkingPrecision]];
  (*Print["traj: "<>ToString[AbsoluteTiming[*)
    (*trajp = traj[wr];*)
  (*][[1]]]];*)
  Print["traj.: "<>ToString[AbsoluteTiming[
  cos\[Chi]r = Cos[\[Chi]r];
  sin\[Chi]r = Sqrt[1-cos\[Chi]r^2];
  rp = p/(1+e*cos\[Chi]r);
  d\[Lambda]d\[Chi]r = Sqrt[(1-e^2)/((1-En0^2) (-p+r3+e r3 cos\[Chi]r) (-p+r4+e r4 cos\[Chi]r))];
  Urp = p*e*sin\[Chi]r/(1+e cos\[Chi]r)^2/d\[Lambda]d\[Chi]r;(* Geodesic radial velocity *)
  \[CapitalDelta]tr = orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][\[Chi]r];
  \[CapitalDelta]\[Phi]r = orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][\[Chi]r];
  \[CapitalDelta]\[Tau]r = orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]r"][\[Chi]r];
  ][[1]]]];
  expr = Exp[I*(\[Omega]*\[CapitalDelta]tr - m*\[CapitalDelta]\[Phi]r + n*\[Chi]r)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
  expr1 = I*(\[Omega]1*\[CapitalDelta]tr + \[Omega]*(-3/(2*Sqrt[K0]))*\[CapitalDelta]\[Tau]r);(* Linear part of the exponential term *)
  Print["Rp: "<>ToString[AbsoluteTiming[
    Rp = R[rp];(* Evaluate the radial functions and their derivatives *)
  ][[1]]]];
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKdr = 2*rp*\[Omega];
  dKd\[Omega]  = (rp^2+a^2);
  d2Kdrd\[Omega] = 2*rp;
  d\[CapitalDelta] = 2*(rp-1);
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  (* In solutions of the radial equation and their r and \[Omega] derivatives *)
  RIn    = Rp["R0"]["In"]["R"]; 
  dRIndr   = Rp["R0"]["In"]["dR"];
  d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  dRInd\[Omega]   = Rp["R1"]["In"]["R"]; (* \[Omega]-derivatives of the radial function *)
  d2RIndrd\[Omega]  = Rp["R1"]["In"]["dR"];
  d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRIn = dRIndr - I*K/\[CapitalDelta]*RIn;
  dDRInd\[Omega] = d2RIndrd\[Omega] - I*K/\[CapitalDelta]*dRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RIn;
  dDRIndr = d2RIndr2 - I*K/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  d2DRIndrd\[Omega] = d3RIndr2d\[Omega] - I*K/\[CapitalDelta]*d2RIndrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRInd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
  dDDRInd\[Omega] = d2DRIndrd\[Omega] - I*K/\[CapitalDelta]*dDRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRIn;
  (* Up solutions of the radial equation and their r and \[Omega] derivatives *)
  RUp    = Rp["R0"]["Up"]["R"]; 
  dRUpdr   = Rp["R0"]["Up"]["dR"];
  d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  dRUpd\[Omega]   = Rp["R1"]["Up"]["R"]; (* \[Omega]-derivative of the radial function *)
  d2RUpdrd\[Omega]  = Rp["R1"]["Up"]["dR"];
  d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRUp = dRUpdr - I*K/\[CapitalDelta]*RUp;
  dDRUpd\[Omega] = d2RUpdrd\[Omega] - I*K/\[CapitalDelta]*dRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RUp;
  dDRUpdr = d2RUpdr2 - I*K/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  d2DRUpdrd\[Omega] = d3RUpdr2d\[Omega] - I*K/\[CapitalDelta]*d2RUpdrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRUpd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
  dDDRUpd\[Omega] = d2DRUpdrd\[Omega] - I*K/\[CapitalDelta]*dDRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRUp;
  (* components of the velocity *)
  vnplus = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
  vnminus = -((rp^2+a^2)*En0 - a*Lz0 - Urp)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
  vmb = (-I*(a*En0 - Lz0));
  (* In solution *)
  {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = FabEq[rp,a,RIn,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
  {dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RIn,DRIn,DDRIn,dRInd\[Omega],dDRInd\[Omega],dDDRInd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
  (* Up solution *)
  {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = FabEq[rp,a,RUp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
  {dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega]}   = dFabd\[Omega]Eq[rp,a,RUp,DRUp,DDRUp,dRUpd\[Omega],dDRUpd\[Omega],dDDRUpd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
  sumPlus0 = Total[d\[Lambda]d\[Chi]r*((vnplus*vnplus*FnnIn + vnplus*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr + (vnminus*vnminus*FnnIn + vnminus*vmb*FnmbIn + vmb*vmb*FmbmbIn)/expr),"Method"->"CompensatedSummation"];
  sumPlus1 = Total[d\[Lambda]d\[Chi]r*((
      (vnplus*GnIn + vmb*GmbIn)/Sqrt[K0] +
      (vnplus*vnplus*FnnIn + vnplus*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr1 +
      (vnplus*vnplus*dFnnInd\[Omega] + vnplus*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
   )*expr + 
   (
      (vnminus*GnIn + vmb*GmbIn)/Sqrt[K0] +
      (vnminus*vnminus*FnnIn + vnminus*vmb*FnmbIn + vmb*vmb*FmbmbIn)*(-expr1) +
      (vnminus*vnminus*dFnnInd\[Omega] + vnminus*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
   )/expr),"Method"->"CompensatedSummation"];
  sumMinus0 = Total[d\[Lambda]d\[Chi]r*((vnplus*vnplus*FnnUp + vnplus*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr + (vnminus*vnminus*FnnUp + vnminus*vmb*FnmbUp + vmb*vmb*FmbmbUp)/expr),"Method"->"CompensatedSummation"];
  sumMinus1 = Total[d\[Lambda]d\[Chi]r*((
      (vnplus*GnUp + vmb*GmbUp)/Sqrt[K0] +
      (vnplus*vnplus*FnnUp + vnplus*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr1 +
      (vnplus*vnplus*dFnnUpd\[Omega] + vnplus*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
   )*expr + 
   (
      (vnminus*GnUp + vmb*GmbUp)/Sqrt[K0] +
      (vnminus*vnminus*FnnUp + vnminus*vmb*FnmbUp + vmb*vmb*FmbmbUp)*(-expr1) +
      (vnminus*vnminus*dFnnUpd\[Omega] + vnminus*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
   )/expr),"Method"->"CompensatedSummation"];
  W = (RIn[[1]]*dRUpdr[[1]] - dRIndr[[1]]*RUp[[1]])/(rp[[1]]^2-2*rp[[1]]+a^2); (* Invariant Wronskian *)
  W1 = (dRInd\[Omega][[1]]*dRUpdr[[1]] + RIn[[1]]*d2RUpdrd\[Omega][[1]] - d2RIndrd\[Omega][[1]]*RUp[[1]] - dRIndr[[1]]*dRUpd\[Omega][[1]])/(rp[[1]]^2-2*rp[[1]]+a^2)*\[Omega]1; (* derivative of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0*\[CapitalOmega]r/(W*stepsr); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0*\[CapitalOmega]r/(W*stepsr);
  CPlus1  = 2*Pi*sumPlus1*\[CapitalOmega]r/(W*stepsr) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CPlus0; (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*sumMinus1*\[CapitalOmega]r/(W*stepsr) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CMinus0;
  \[ScriptCapitalF]En\[ScriptCapitalI] = Abs[CPlus0]^2/(4Pi*\[Omega]^2); (* Fluxes and their linear parts *)
  \[ScriptCapitalF]En\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI] = Abs[CPlus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH] = \[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]En\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]En\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2);
  \[ScriptCapitalF]Lz\[ScriptCapitalI]1 = (2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  \[ScriptCapitalF]Lz\[ScriptCapitalH]1 = \[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3);
  <|
    "l"->l, "m"->m, "k"->0, "n"->n,
    "\[Omega]"->\[Omega], "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"-><|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"-><|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha], "\[Alpha]Correction"->\[Alpha]1,
    "S"->S,
    "Fluxes"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI], "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]|>
    |>,
    "FluxesCorrection"-><|
      "Energy"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]En\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]En\[ScriptCapitalH]1|>,
      "AngularMomentum"-><|"\[ScriptCapitalI]"->\[ScriptCapitalF]Lz\[ScriptCapitalI]1, "\[ScriptCapitalH]"->\[ScriptCapitalF]Lz\[ScriptCapitalH]1|>
    |>,
    "stepsr"->stepsr
  |>
]


(* ::Subsection::Closed:: *)
(*Generic*)


(* ::Subsubsection::Closed:: *)
(*Numerical trajectory*)


TeukolskySpinMode[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbit_]:=Module[{a,p,e,\[ScriptCapitalI],s,En0,Lz0,Kc0,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalGamma],\[Omega],SWSH,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Theta]list,
    R,\[Lambda],W,sumPlus,sumMinus,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Urp,Uzp,udtp,ud\[Phi]p,Kp,\[CapitalDelta],dK,d\[CapitalDelta],\[CapitalSigma],\[Zeta],\[Zeta]bar,fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,
    dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,C0nn,C0nmb,C0mbmb,
    Crnn,Crnmb,Crmbmb,C\[Theta]nn,C\[Theta]nmb,C\[Theta]mbmb,Amnn,Amnmb,Ammbmb,Adnn,Adnmb,Admbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,
    A0,A1,A2,B1,B2,B3,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,RUprp,dRUprp,ddRUprp,dddRUprp,CPlus,CMinus,\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],expr,exp\[Theta],exp1,stepsr,steps\[Theta],h1,h2,h3,
    r0,z0,Ur0,Uz0,r0p,z0p,\[CapitalDelta]0,d\[CapitalDelta]0,K0,dK0,r1p,z1p,sin\[Theta]0p,\[Zeta]0,\[Zeta]0bar,Ur1p,Uz1p,correction,correctionp},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Kc0 (r^2-3 a^2 z^2)))/(Sqrt[Kc0] (r^2+a^2 z^2)^3);(* Projections of s^\[Mu]\[Nu] to \[Xi]_\[Mu];\[Nu] in Eq. (21) in [2303.16798] *)
  h2[r_,z_] := 1/(Sqrt[Kc0] (r^2+a^2 z^2)^3) (-En0 Lz0 r^6+a^4 En0 Lz0 r^2 z^4+a^2 En0 Lz0 r^4 (-2+z^2)-a^6 En0 Lz0 z^4 (-2+z^2)+
               a^7 En0^2 z^4 (-1+z^2)+a r^3 (Lz0^2 r+Kc0 (-1+z^2)+Kc0 r (-1+2 z^2)+r^3 (z^2-En0^2 (-1+z^2)))+a^3 (-En0^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc0 r z^2-3 Kc0 (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc0-Lz0^2+r ( (-1+z^2)+r (2-z^2+En0^2 (-1+z^2)))));
  h3[r_,z_]:=-((2 a r z)/(Sqrt[Kc0] (r^2+a^2 z^2)^2));
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  \[ScriptCapitalI] = orbit["\[ScriptCapitalI]"];
  s = orbit["s"];
  En0 = orbit["Ehat"];(* Geodesic constants of motion *)
  Lz0 = orbit["Lzhat"];
  Kc0 = orbit["Khat"];
  En1 = orbit["ES"];(* Linear corrections to the constants of motion *)
  Lz1 = orbit["LS"];
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbit["BLFrequencies"]; (* BL frequencies *)
  correction = orbit["OrbitCorrection"];(* function containing corrections to the trajectory *)
  r0 = orbit["TrajectoryGeo"][[2]];(* Geodesic coordinates r and z=cos(\[Theta]) *)
  z0[wz_]:=Cos[orbit["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbit["MinoFrequencies"][[1]];(* Average rate of change of BL time in Mino time *)
  \[Omega] = m*\[CapitalOmega]\[Phi]+n*\[CapitalOmega]r+k*\[CapitalOmega]\[Theta]; (* Frequency *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue *)
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega],
        Method->{"NumericalIntegration","Domain"->{"In"->{p/(1+e),p/(1-e)},"Up"->{p/(1+e),p/(1-e)}}}];
  \[Lambda] = R["In"]["Eigenvalue"];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);  (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus = 0; (* Results of the integration are stored in these variables *)
  sumMinus = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in w\[Theta]"];
  \[Theta]list = {};(* List for functions of \[Theta] *)
  correctionp = Table[correction[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]],N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]]],{ir,1,stepsr/2},{i\[Theta],1,steps\[Theta]/2}];(* Corrections to the trajectory *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]];
    r0p = r0[wr];(* Geodesic r *)
    Ur0 = {1,-1,1,-1}*Sqrt[((r0p^2+a^2)*En0-a*Lz0)^2-(r0p^2-2r0p+a^2)*(r0p^2+Kc0)];(* Geodesic radial velocity at each quadrant *)
    expr = Exp[I*(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr]+2Pi*n*(ir-1/2)/stepsr)];
    \[CapitalDelta]0 = r0p^2-2r0p+a^2;
    K0 = (r0p^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta]0 = 2*(r0p-1);
    dK0 = 2r0p*\[Omega];
    V = -(K0^2+4I*(r0p-1)*K0)/\[CapitalDelta]0+8I*\[Omega]*r0p+\[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2K0*dK0+4I*K0+4I*(r0p-1)*dK0)*\[CapitalDelta]0-(K0^2+4I*(r0p-1)*K0)*d\[CapitalDelta]0)/\[CapitalDelta]0^2+8I*\[Omega]; (* derivative of the potential in radial Teukolsky equation *)
    RInrp = R["In"][r0p]; (* radial function *)
    dRInrp = R["In"]'[r0p];
    ddRInrp = (V*RInrp+d\[CapitalDelta]0*dRInrp)/\[CapitalDelta]0;  (* second derivative of radial function from Teukolsky equation *)
    dddRInrp = 1/\[CapitalDelta]0 (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp = R["Up"][r0p];
    dRUprp = R["Up"]'[r0p];
    ddRUprp = (V*RUprp+d\[CapitalDelta]0*dRUprp)/\[CapitalDelta]0;  
    dddRUprp = 1/\[CapitalDelta]0 (dV*RUprp+(V+2)*dRUprp);
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++,(* integral over w_\[Theta] *)
      w\[Theta] = N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]];
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        z0p = z0[w\[Theta]];(* Geodesic z=cos(\[Theta]) *)
        Uz0 = {1,1,-1,-1}*(-1)*Sqrt[-((1-z0p^2)*a*En0-Lz0)^2+(1-z0p^2)*(Kc0-a^2*z0p^2)];(* Geodesic polar velocity at each quadrant *)
        exp\[Theta] = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]0p = Sqrt[1-z0p^2];
        S = SWSH[ArcCos[z0p],0];  (*  Angular function S(\[Theta])  *)
        dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[z0p]);
        d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-z0p^2)-(m-2*z0p)^2/(1-z0p^2)+4a*\[Omega]*z0p-2+2*m*a*\[Omega]+\[Lambda])*S-z0p/sin\[Theta]0p*dSd\[Theta];(* second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]0p^3)2 (-2+m z0p+a z0p (-1+z0p^2) \[Omega]) (m+a \[Omega]-z0p (2+a z0p \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-z0p^2)-(m-2*z0p)^2/(1-z0p^2)+4a*\[Omega]*z0p-2+2*m*a*\[Omega]+\[Lambda]-1/(1-z0p^2))*dSd\[Theta]-z0p/sin\[Theta]0p*d2Sd\[Theta]2; (*third derivative of S from Teukolsky equation *)
        L2S = dSd\[Theta]-(S (m-2 z0p+a (-1+z0p^2) \[Omega]))/sin\[Theta]0p;
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+z0p^2) (dSd\[Theta] sin\[Theta]0p (m-2 z0p+a (-1+z0p^2) \[Omega])+S (2-m z0p+a z0p (-1+z0p^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 z0p-2 a (-1+z0p^2) \[Omega]))/sin\[Theta]0p+S (-2-(m (m-2 z0p))/(-1+z0p^2)-2 a (m-2 z0p) \[Omega]-a^2 (-1+z0p^2) \[Omega]^2);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+z0p^2) dSd\[Theta] (5-m^2-2 z0p^2-2 a (m-3 z0p) (-1+z0p^2) \[Omega]-a^2 (-1+z0p^2)^2 \[Omega]^2)+
                   1/sin\[Theta]0p^3 (d2Sd\[Theta]2 (-1+z0p^2) (2 m-3 z0p+2 a (-1+z0p^2) \[Omega])+2 S (-m^2 z0p+m (1+z0p^2)+a (-1+z0p^2)^2 \[Omega] (-2+a z0p \[Omega])));
        AppendTo[\[Theta]list,{z0p,Uz0,exp\[Theta],sin\[Theta]0p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}];,
        {z0p,Uz0,exp\[Theta],sin\[Theta]0p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]} = \[Theta]list[[i\[Theta]]];
      ];
      \[Zeta]0 = r0p-I*a*z0p;
      \[Zeta]0bar = r0p+I*a*z0p;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta]0,\[Zeta]0bar,sin\[Theta]0p,\[CapitalDelta]0,d\[CapitalDelta]0,K0,dK0,S,L2S,L1L2S,a];(* Functions f_ab^(i) and their derivatives *)
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta]0,\[Zeta]0bar,sin\[Theta]0p,\[CapitalDelta]0,d\[CapitalDelta]0,K0,dK0,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta]0,\[Zeta]0bar,sin\[Theta]0p,\[CapitalDelta]0,d\[CapitalDelta]0,K0,dK0,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      r1p = {correctionp[[ir,i\[Theta]]]["rS"], correctionp[[ir,-i\[Theta]]]["rS"], correctionp[[ir,-i\[Theta]]]["rS"],correctionp[[ir,i\[Theta]]]["rS"]};(* Corrections to the coordiantes and four-velocities at each quadrant *)
      z1p = {correctionp[[ir,i\[Theta]]]["zS"],-correctionp[[ir,-i\[Theta]]]["zS"],-correctionp[[ir,-i\[Theta]]]["zS"],correctionp[[ir,i\[Theta]]]["zS"]};
      Ur1p = {correctionp[[ir,i\[Theta]]]["UrS"],-correctionp[[ir,-i\[Theta]]]["UrS"], correctionp[[ir,-i\[Theta]]]["UrS"],-correctionp[[ir,i\[Theta]]]["UrS"]};
      Uz1p = {correctionp[[ir,i\[Theta]]]["UzS"], correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,i\[Theta]]]["UzS"]};
      rp = r0p + s*r1p;
      zp = z0p + s*z1p;
      sin\[Theta]p = Sqrt[1-zp^2];
      Urp = Ur0 + s*Ur1p;
      Uzp = Uz0 + s*Uz1p;
      udtp = -(En0+s*En1) + s*h1[r0p,z0p];(* covariant t and \[Phi] components of the four-velocity *)
      ud\[Phi]p = (Lz0+s*Lz1) + s*(h2[r0p,z0p] + h3[r0p,z0p]*Ur0*Uz0);
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      Kp = (rp^2+a^2)*\[Omega]-a*m;
      \[CapitalDelta] = rp^2-2*rp+a^2;
      dK = 2*rp*\[Omega];
      d\[CapitalDelta] = 2*(rp-1);
      \[CapitalSigma] = rp^2+a^2*zp^2;
      exp1 = {Exp[ I*s*(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])],
              Exp[-I*s*(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])],
              Exp[ I*s*(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])],
              Exp[-I*s*(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])]};(* Correction to the exponential term *)
      vn = ((rp^2+a^2)*udtp+a*ud\[Phi]p-Urp)/(2*\[CapitalSigma]);(* Components of the four-velocity in the Kinnersley tetrad *)
      vl = ((rp^2+a^2)*udtp+a*ud\[Phi]p+Urp)/(\[CapitalDelta]);
      vm = (-I*(a*sin\[Theta]p^2*udtp+ud\[Phi]p)+Uzp)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      Sln = s*(-((rp (Kc0-a^2 zp^2))/(Sqrt[Kc0] \[CapitalSigma])));(* Spin tensor in Kinnersley tetrad *)
      Snm = s*(\[Zeta]/Sqrt[Kc0])*vm*vn;
      Snmb = Conjugate[Snm];
      Slmb = s*(-(\[Zeta]/Sqrt[Kc0]))*vl*vmb;
      Smmb = s*((I a zp (Kc0+rp^2))/(Sqrt[Kc0] \[CapitalSigma]));
      Amnn   = vn^2;
      Amnmb  = vn*vmb;
      Ammbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);(* term in the comment is indentically zero *)
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Adnn   = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Admbmb = (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Adnmb  = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)-Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*Kp/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*Kp*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      C0nn   = (Amnn+Adnn+St\[Phi]n*vn);
      Crnn   = (Amnn*s*r1p+Srn*vn);
      C\[Theta]nn   = (-Amnn*s*z1p/sin\[Theta]0p+S\[Theta]n*vn);
      C0nmb  = (Amnmb+Adnmb+(St\[Phi]n*vmb+St\[Phi]mb*vn)/2);
      Crnmb  = (Amnmb*s*r1p+(Srn*vmb+Srmb*vn)/2);
      C\[Theta]nmb  = (-Amnmb*s*z1p/sin\[Theta]0p+(S\[Theta]n*vmb+S\[Theta]mb*vn)/2);
      C0mbmb = (Ammbmb+Admbmb+St\[Phi]mb*vmb);
      Crmbmb = (Ammbmb*s*r1p+Srmb*vmb);
      C\[Theta]mbmb = (-Ammbmb*s*z1p/sin\[Theta]0p+S\[Theta]mb*vmb);
      A0 = C0nn*fnn0 + Crnn*dfnn0dr + C\[Theta]nn*dfnn0d\[Theta] + C0nmb*fnmb0 + Crnmb*dfnmb0dr + C\[Theta]nmb*dfnmb0d\[Theta] + C0mbmb*fmbmb0 + Crmbmb*dfmbmb0dr + C\[Theta]mbmb*dfmbmb0d\[Theta];
      A1 = C0nmb*fnmb1 + Crnmb*dfnmb1dr + C\[Theta]nmb*dfnmb1d\[Theta] + C0mbmb*fmbmb1 + Crmbmb*dfmbmb1dr + C\[Theta]mbmb*dfmbmb1d\[Theta];
      A2 = C0mbmb*fmbmb2 + Crmbmb*dfmbmb2dr + C\[Theta]mbmb*dfmbmb2d\[Theta];
      B1 = - Srn*vn*fnn0 - (Srn*vmb+Srmb*vn)/2*fnmb0 - Srmb*vmb*fmbmb0;
      B2 = - (Srn*vmb+Srmb*vn)/2*fnmb1 - Srmb*vmb*fmbmb1;
      B3 = - Srmb*vmb*fmbmb2;
      sumPlus += Total[\[CapitalSigma]*((A0*(RInrp+s*r1p*dRInrp)-(A1+B1)*(dRInrp+s*r1p*ddRInrp)+(A2+B2)*(ddRInrp+s*r1p*dddRInrp)-B3*dddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}*exp1)]; (* Total over all quadrants *)
      sumMinus += Total[\[CapitalSigma]*((A0*(RUprp+s*r1p*dRUprp)-(A1+B1)*(dRUprp+s*r1p*ddRUprp)+(A2+B2)*(ddRUprp+s*r1p*dddRUprp)-B3*dddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}*exp1)];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/(r0p^2-2*r0p+a^2); (* Invariant Wronskian *)
  CPlus = 2Pi/\[CapitalGamma]/W/stepsr/steps\[Theta]*sumPlus; (* Amplitudes *)
  CMinus = 2Pi/\[CapitalGamma]/W/stepsr/steps\[Theta]*sumMinus;
  Association[
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus,
      "\[ScriptCapitalH]"->CMinus
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  ] (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


TeukolskySpinModeFromCorrection[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbitCorrection_,s_]:=Module[{h1,h2,h3,a,p,e,\[ScriptCapitalI],En0,Lz0,Kc0,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],correction,r,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,
    sumPlus,sumMinus,stepsr,steps\[Theta],\[Theta]list,correctionp,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,
    RUprp,dRUprp,ddRUprp,dddRUprp,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,
    dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,
    Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,r1p,z1p,Ur1p,Uz1p,\[CapitalSigma]1,exp1,vn1,vmb1,
    Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus,CMinus},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc0 (r^2-3 a^2 z^2)))/(Sqrt[Kc0] (r^2+a^2 z^2)^3);(* Projections of s^\[Mu]\[Nu] to \[Xi]_\[Mu];\[Nu] in Eq. (21) in [2303.16798] *)
  h2[r_,z_] := 1/(Sqrt[Kc0] (r^2+a^2 z^2)^3) (-En0 Lz0 r^6+a^4 En0 Lz0 r^2 z^4+a^2 En0 Lz0 r^4 (-2+z^2)-a^6 En0 Lz0 z^4 (-2+z^2)+
               a^7 En0^2 z^4 (-1+z^2)+a r^3 (Lz0^2 r+Kc0 (-1+z^2)+Kc0 r (-1+2 z^2)+r^3 (z^2-En0^2 (-1+z^2)))+a^3 (-En0^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc0 r z^2-3 Kc0 (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc0-Lz0^2+r ( (-1+z^2)+r (2-z^2+En0^2 (-1+z^2)))));
  h3[r_,z_] := -((2 a r z)/(Sqrt[Kc0] (r^2+a^2 z^2)^2));
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  \[ScriptCapitalI] = orbitCorrection["\[ScriptCapitalI]"];
  En0 = orbitCorrection["Ehat"]; (* Geodesic constants of motion *)
  Lz0 = orbitCorrection["Lzhat"];
  Kc0 = orbitCorrection["Khat"];
  En1 = orbitCorrection["ES"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["LS"];
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]+s*orbitCorrection["BLFrequenciesCorrection"]; (* BL frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function computing corrections to the trajectory *)
  r = orbitCorrection["TrajectoryGeo"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Frequency of mode *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue *)
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega],Method->{"NumericalIntegration",
        "Domain"->{"In"->{p/(1+e),p/(1-e)},"Up"->{p/(1+e),p/(1-e)}}}];
  \[Lambda] = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalEigenvalue[-2,l,m,a*\[Omega]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);  (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (* Outer horizon *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus  = 0; (* Results of the integration are stored in these variables *)
  sumMinus = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in w\[Theta]"];
  \[Theta]list = {};(* List for functions of \[Theta] *)
  correctionp = Table[correction[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]],N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]]],{ir,1,stepsr/2},{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]];
    rp = r[wr];
    Ur = {1,-1,1,-1}*If[wr<Pi,1,-1]*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+Kc0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2rp*\[Omega];
    V  = -(K^2+4I*(rp-1)*K)/\[CapitalDelta]+8I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation *)
    RInrp    = R["In"][rp]; (* radial function *)
    dRInrp   = R["In"]'[rp];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R["Up"][rp];
    dRUprp   = R["Up"]'[rp];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++,(* integration over w\[Theta] *)
      w\[Theta] = N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]];
      If[ir==1,(* functions of only \[Theta] saved to a list in the first step *)
        zp = z[w\[Theta]];
        Uz = {1,1,-1,-1}*If[w\[Theta]<Pi,-1,1]*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(Kc0-a^2*zp^2)];(* Polar geodesic velocity *)
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2];
        S = SWSH[ArcCos[zp],0];  (*  Polar function S(\[Theta](z))  *)
        dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]);
        d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta];(* second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}],
        {zp,Uz,exp\[Theta],sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      vn = -((rp^2+a^2)*En0 - a*Lz0 + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En0 - a*Lz0 - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En0 - Lz0) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      Sln  = (-((rp (Kc0-a^2 zp^2))/(Sqrt[Kc0] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
      Snm  = (\[Zeta]/Sqrt[Kc0])*vm*vn;
      Snmb = Conjugate[Snm];
      Slmb = (-(\[Zeta]/Sqrt[Kc0]))*vl*vmb;
      Smmb = ((I a zp (Kc0+rp^2))/(Sqrt[Kc0] \[CapitalSigma]));
      Amnn   = vn^2; (* Monopole terms *)
      Amnmb  = vn*vmb;
      Ammbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients (with opposite sign than in Teukolsky paper) *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Adnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]); (* Dipole terms *)
      Admbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Adnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm); (* Functions from Eq. (B3) in [2303.16798] *)
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      r1p = {correctionp[[ir,i\[Theta]]]["rS"], correctionp[[ir,-i\[Theta]]]["rS"], correctionp[[ir,-i\[Theta]]]["rS"],correctionp[[ir,i\[Theta]]]["rS"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
      z1p = {correctionp[[ir,i\[Theta]]]["zS"],-correctionp[[ir,-i\[Theta]]]["zS"],-correctionp[[ir,-i\[Theta]]]["zS"],correctionp[[ir,i\[Theta]]]["zS"]};
      Ur1p = {correctionp[[ir,i\[Theta]]]["UrS"],-correctionp[[ir,-i\[Theta]]]["UrS"], correctionp[[ir,-i\[Theta]]]["UrS"],-correctionp[[ir,i\[Theta]]]["UrS"]};
      Uz1p = {correctionp[[ir,i\[Theta]]]["UzS"], correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,i\[Theta]]]["UzS"]};
      \[CapitalSigma]1 = 2*(rp*r1p+a^2*zp*z1p); (* Linear correction to \[CapitalSigma] *)
      exp1 = {I*(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"]),-I*(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"]),
              I*(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"]),-I*(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])};(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      vn1  = -( ((2*rp*r1p)*En0) - ((rp^2+a^2)*En0 - a*Lz0 + Ur)/(\[CapitalSigma])*\[CapitalSigma]1 + 
               ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz) + Ur1p))/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = ( (-I*(-2*a*zp*z1p*En0)) - (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uz)*(-zp*z1p/sin\[Theta]p^2 + (r1p-I*a*z1p)/\[Zeta]) + 
               (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz)) + Uz1p))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
      Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Adnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + r1p*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - z1p/sin\[Theta]p*vn)*vn;
      Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Adnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + r1p*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - z1p/sin\[Theta]p*vn*vmb);
      Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Admbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + r1p*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - z1p/sin\[Theta]p*vmb)*vmb;
      sumPlus  += Total[s*\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RInrp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRInrp + (Ammbmb*fmbmb2)*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})];  (* Total over all quadrants *)
      sumMinus += Total[s*\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RUprp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRUprp + (Ammbmb*fmbmb2)*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  CPlus  = 2*Pi*(1-s*\[CapitalGamma]1/\[CapitalGamma])*sumPlus/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Amplitudes *)
  CMinus = 2*Pi*(1-s*\[CapitalGamma]1/\[CapitalGamma])*sumMinus/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus,
      "\[ScriptCapitalH]"->CMinus
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


dJzdK[a_,s_,{Eso_,Lso_,Kso_}]:=Module[{z1,z2,kzsq,dJz0dKso,dJz1dKso},
z2 = Sqrt[1/(2 (a^2-a^2 Eso^2)) (a^2-2 a^2 Eso^2+Kso+2 a Eso Lso-Sqrt[a^4-2 a^2 Kso+Kso^2-4 a^3 Eso Lso+4 a Eso Kso Lso+4 a^2 Lso^2])];
z1 = Sqrt[1/(2 (a^2-a^2 Eso^2)) (a^2-2 a^2 Eso^2+Kso+2 a Eso Lso+Sqrt[a^4-2 a^2 Kso+Kso^2-4 a^3 Eso Lso+4 a Eso Kso Lso+4 a^2 Lso^2])];
kzsq=z2^2/z1^2;
dJz0dKso=EllipticK[kzsq]/(a Sqrt[1-Eso^2] \[Pi] z1);
dJz1dKso=-(1/(a Sqrt[1-Eso^2] Sqrt[Kso] \[Pi] z1)) (Eso EllipticK[kzsq]+((a^2 Eso-Eso Kso-a Lso) (a^2 z1^2 EllipticE[kzsq]+(Kso-a^2 z1^2) EllipticK[kzsq]))/((Kso-a^2 z1^2) (Kso-a^2 z2^2))-(Kso (-(((1-z1^2) (Lso-a Eso (1-z1^2)) EllipticE[kzsq])/((1-kzsq) z1^2 (Kso-a^2 z1^2)))+((1-z2^2) (Lso-a Eso (1-z2^2)) (-(EllipticE[kzsq]/(1-kzsq))+EllipticK[kzsq]))/(z2^2 (Kso-a^2 z2^2))))/(a (1-Eso^2) (-z1^2+z2^2)));
dJz0dKso + s*dJz1dKso
];

dJzdEn[a_,s_,{Eso_,Lso_,Kso_}]:=Module[{z1,z2,kzsq,dJz0dEso,dJz1dEso},
z2 = Sqrt[1/(2 (a^2-a^2 Eso^2)) (a^2-2 a^2 Eso^2+Kso+2 a Eso Lso-Sqrt[a^4-2 a^2 Kso+Kso^2-4 a^3 Eso Lso+4 a Eso Kso Lso+4 a^2 Lso^2])];
z1 = Sqrt[1/(2 (a^2-a^2 Eso^2)) (a^2-2 a^2 Eso^2+Kso+2 a Eso Lso+Sqrt[a^4-2 a^2 Kso+Kso^2-4 a^3 Eso Lso+4 a Eso Kso Lso+4 a^2 Lso^2])];
kzsq=z2^2/z1^2;
dJz0dEso=1/(Sqrt[1-Eso^2] \[Pi] z1) 2 (-a Eso z1^2 EllipticE[kzsq]+(Lso-a Eso (1-z1^2)) EllipticK[kzsq]);
dJz1dEso=(2 Sqrt[Kso] ((z1^2+z2^2 (1+z1^4+z1^2 (-4+z2^2))) EllipticE[kzsq]+(-z1+z2) (z1+z2) (1+(-2+z1^2) z2^2) EllipticK[kzsq]))/(a (1-Eso^2)^(3/2) \[Pi] z1 z2^2 (-z1^2+z2^2)^2);
dJz0dEso + s*dJz1dEso
];

dJzdLz[a_,s_,{Eso_,Lso_,Kso_}]:=Module[{z1,z2,kzsq,dJz0dLso,dJz1dLso},
z2 = Sqrt[1/(2 (a^2-a^2 Eso^2)) (a^2-2 a^2 Eso^2+Kso+2 a Eso Lso-Sqrt[a^4-2 a^2 Kso+Kso^2-4 a^3 Eso Lso+4 a Eso Kso Lso+4 a^2 Lso^2])];
z1 = Sqrt[1/(2 (a^2-a^2 Eso^2)) (a^2-2 a^2 Eso^2+Kso+2 a Eso Lso+Sqrt[a^4-2 a^2 Kso+Kso^2-4 a^3 Eso Lso+4 a Eso Kso Lso+4 a^2 Lso^2])];
kzsq=z2^2/z1^2;
dJz0dLso=-((2 (-a Eso EllipticK[kzsq]+Lso EllipticPi[z2^2,kzsq]))/(a Sqrt[1-Eso^2] \[Pi] z1));
dJz1dLso=(2 Sqrt[Kso] ((-z1^2+(-1+2 z1^2) z2^2) EllipticE[kzsq]+(-z1+z2) (z1+z2) (-1+z2^2) EllipticK[kzsq]))/(a^2 (1-Eso^2)^(3/2) \[Pi] z1 z2^2 (-z1^2+z2^2)^2);
dJz0dLso + s*dJz1dLso
];

KCorrection[a_,p_,e_,x_,{Eng_,Lzg_,Kg_}]:=Module[{IntR\[Xi]t\[Phi],IntZ\[Xi]t\[Phi],az12Int\[CapitalSigma],dRdEn,dRdLz,dRdK,dZdEn,dZdLz,dZdK,r1,r2,r3,r4,z1,z2,kz,kr,v,M,ES,LzS,KS,r1S,r2S,z1S,r},
IntR\[Xi]t\[Phi][r_,En_,Lz_,K_,z1_,kz_]:=-2(a^2 En-a Lz+En r^2)/Sqrt[K]*((r*((r-1)*(K+r^2)+r*(r^2-2r+a^2))-2*r^2*En*((r^2+a^2)*En-a*Lz))*EllipticPi[-a^2*z1^2/r^2,kz]/EllipticK[kz]/r^2+(En*((r^2+a^2)*En-a*Lz)-r^2+r));
dRdEn[r_,En_,Lz_,K_]:=2 (a^2+r^2) (-a Lz+En (a^2+r^2));
dRdLz[r_,En_,Lz_,K_]:=-2 a (-a Lz+En (a^2+r^2));
dRdK[r_,En_,Lz_,K_]:=-a^2-(-2+r) r;
dZdEn[z_,En_,Lz_,K_]:=-2 a (-1+z^2) (Lz+a En (-1+z^2));
dZdLz[z_,En_,Lz_,K_]:=-2 (Lz+a En (-1+z^2));
dZdK[z_,En_,Lz_,K_]:=1-z^2;
{r1,r2,r3,r4}=Reverse[Sort[NSolveValues[((r^2+a^2)*Eng-a*Lzg)^2-(r^2-2r+a^2)*(Kg+r^2)==0,r]]];
z1=Sqrt[1-x^2];
z2=Sqrt[a^2*(1-Eng^2)+Lzg^2/x^2];
kz=a^2*(1-Eng^2)*z1^2/z2^2;
kr=(r1-r2)*(r3-r4)/((r1-r3)*(r2-r4));
az12Int\[CapitalSigma]=((a*z1^2)/(r3^2+a^2 z1^2)+2*Re[((I (-r2+r3)z1)/(2 (r2+I a z1) (r3+I a z1)))*EllipticPi[((r1-r2) (r3+I a z1))/((r1-r3) (r2+I a z1)),kr]]/EllipticK[kr]);
IntZ\[Xi]t\[Phi]=-2*(Lzg+a Eng (-1+z1^2))/Sqrt[Kg]*((Eng*((1-z1^2)*a*Eng-Lzg)+a z1^2)+((Kg+a^2 (1-2 z1^2))-2*a*Eng*((1-z1^2)*a*Eng-Lzg))*az12Int\[CapitalSigma]);
v={IntR\[Xi]t\[Phi][r1,Eng,Lzg,Kg,z1,kz],IntR\[Xi]t\[Phi][r2,Eng,Lzg,Kg,z1,kz],IntZ\[Xi]t\[Phi]};
M={{dRdEn[r1,Eng,Lzg,Kg],dRdLz[r1,Eng,Lzg,Kg],dRdK[r1,Eng,Lzg,Kg]},{dRdEn[r2,Eng,Lzg,Kg],dRdLz[r2,Eng,Lzg,Kg],dRdK[r2,Eng,Lzg,Kg]},{dZdEn[z1,Eng,Lzg,Kg],dZdLz[z1,Eng,Lzg,Kg],dZdK[z1,Eng,Lzg,Kg]}};
{ES,LzS,KS}=LinearSolve[M,-v];

KS+2*((Kg+a^2)*Eng-a*Lzg)/Sqrt[Kg]
]


TeukolskySpinModeFromCorrectionNew[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbitCorrection_,s_]:=Module[{h1,h2,h3,a,p,e,\[ScriptCapitalI],En0,Lz0,Kc0,En1,Lz1,Kc1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],correction,r,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,
    sumPlus,sumMinus,stepsr,steps\[Theta],\[Theta]list,correctionp,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,
    RUprp,dRUprp,ddRUprp,dddRUprp,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,
    dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,
    Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,r1p,z1p,Ur1p,Uz1p,\[CapitalSigma]1,exp1,vn1,vmb1,
    Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus,CMinus,FnnIn,FnmbIn,FmbmbIn,
    dFnnIndr,dFnmbIndr,dFmbmbIndr,dFnnInd\[Theta],dFnmbInd\[Theta],dFmbmbInd\[Theta],FnnUp,FnmbUp,FmbmbUp,
    dFnnUpdr,dFnmbUpdr,dFmbmbUpdr,dFnnUpd\[Theta],dFnmbUpd\[Theta],dFmbmbUpd\[Theta]},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc0 (r^2-3 a^2 z^2)))/(Sqrt[Kc0] (r^2+a^2 z^2)^3);(* Projections of s^\[Mu]\[Nu] to \[Xi]_\[Mu];\[Nu] in Eq. (21) in [2303.16798] *)
  h2[r_,z_] := 1/(Sqrt[Kc0] (r^2+a^2 z^2)^3) (-En0 Lz0 r^6+a^4 En0 Lz0 r^2 z^4+a^2 En0 Lz0 r^4 (-2+z^2)-a^6 En0 Lz0 z^4 (-2+z^2)+
               a^7 En0^2 z^4 (-1+z^2)+a r^3 (Lz0^2 r+Kc0 (-1+z^2)+Kc0 r (-1+2 z^2)+r^3 (z^2-En0^2 (-1+z^2)))+a^3 (-En0^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc0 r z^2-3 Kc0 (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc0-Lz0^2+r ( (-1+z^2)+r (2-z^2+En0^2 (-1+z^2)))));
  h3[r_,z_] := -((2 a r z)/(Sqrt[Kc0] (r^2+a^2 z^2)^2));
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  \[ScriptCapitalI] = orbitCorrection["\[ScriptCapitalI]"];
  En0 = orbitCorrection["Ehat"]; (* Geodesic constants of motion *)
  Lz0 = orbitCorrection["Lzhat"];
  Kc0 = orbitCorrection["Khat"];
  En1 = orbitCorrection["ES"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["LS"];
  Kc1 = KCorrection[a,p,e,Cos[\[ScriptCapitalI]],{En0,Lz0,Kc0}];
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]+s*orbitCorrection["BLFrequenciesCorrection"]; (* BL frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function computing corrections to the trajectory *)
  r = orbitCorrection["TrajectoryGeo"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Frequency of mode *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue *)
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega],Method->{"NumericalIntegration",
        "Domain"->{"In"->{p/(1+e),p/(1-e)},"Up"->{p/(1+e),p/(1-e)}}}];
  \[Lambda] = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalEigenvalue[-2,l,m,a*\[Omega]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);  (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (* Outer horizon *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus  = 0; (* Results of the integration are stored in these variables *)
  sumMinus = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in w\[Theta]"];
  \[Theta]list = {};(* List for functions of \[Theta] *)
  correctionp = Table[correction[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]],N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]]],{ir,1,stepsr/2},{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]];
    rp = r[wr];
    Ur = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+Kc0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2rp*\[Omega];
    V  = -(K^2+4I*(rp-1)*K)/\[CapitalDelta]+8I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation *)
    RInrp    = R["In"][rp]; (* radial function *)
    dRInrp   = R["In"]'[rp];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R["Up"][rp];
    dRUprp   = R["Up"]'[rp];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++,(* integration over w\[Theta] *)
      w\[Theta] = N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]];
      If[ir==1,(* functions of only \[Theta] saved to a list in the first step *)
        zp = z[w\[Theta]];
        Uz = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(Kc0-a^2*zp^2)];(* Polar geodesic velocity *)
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2];
        S = SWSH[ArcCos[zp],0];  (*  Polar function S(\[Theta](z))  *)
        dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]);
        d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta];(* second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}],
        {zp,Uz,exp\[Theta],sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      vn = -((rp^2+a^2)*En0 - a*Lz0 + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En0 - a*Lz0 - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En0 - Lz0) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      Sln  = (-((rp (Kc0-a^2 zp^2))/(Sqrt[Kc0] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
      Snm  = (\[Zeta]/Sqrt[Kc0])*vm*vn;
      Snmb = Conjugate[Snm];
      Slmb = (-(\[Zeta]/Sqrt[Kc0]))*vl*vmb;
      Smmb = ((I a zp (Kc0+rp^2))/(Sqrt[Kc0] \[CapitalSigma]));
      rho = 1/\[Zeta]; (* Spin coefficients (with opposite sign than in Teukolsky paper) *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Adnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]); (* Dipole terms *)
      Admbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Adnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm); (* Functions from Eq. (B3) in [2303.16798] *)
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      r1p = {correctionp[[ir,i\[Theta]]]["rS"], correctionp[[ir,-i\[Theta]]]["rS"], correctionp[[ir,-i\[Theta]]]["rS"],correctionp[[ir,i\[Theta]]]["rS"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
      z1p = {correctionp[[ir,i\[Theta]]]["zS"],-correctionp[[ir,-i\[Theta]]]["zS"],-correctionp[[ir,-i\[Theta]]]["zS"],correctionp[[ir,i\[Theta]]]["zS"]};
      Ur1p = {correctionp[[ir,i\[Theta]]]["UrS"],-correctionp[[ir,-i\[Theta]]]["UrS"], correctionp[[ir,-i\[Theta]]]["UrS"],-correctionp[[ir,i\[Theta]]]["UrS"]};
      Uz1p = {correctionp[[ir,i\[Theta]]]["UzS"], correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,i\[Theta]]]["UzS"]};
      \[CapitalSigma]1 = 2*(rp*r1p+a^2*zp*z1p); (* Linear correction to \[CapitalSigma] *)
      exp1 = {I*(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"]),-I*(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"]),
              I*(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"]),-I*(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])};(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      vn1  = -( ((2*rp*r1p)*En0) - ((rp^2+a^2)*En0 - a*Lz0 + Ur)/(\[CapitalSigma])*\[CapitalSigma]1 + 
               ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz) + Ur1p))/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = ( (-I*(-2*a*zp*z1p*En0)) - (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uz)*(-zp*z1p/sin\[Theta]p^2 + (r1p-I*a*z1p)/\[Zeta]) + 
               (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz)) + Uz1p))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
      Amnn   = (vn + s*vn1)^2; (* Monopole terms *)
      Amnmb  = (vn + s*vn1)*(vmb + s*vmb1);
      Ammbmb = (vmb + s*vmb1)^2;
      Ann0S  = Adnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + r1p*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - z1p/sin\[Theta]p*vn)*vn;
      Anmb0S  = (Adnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + r1p*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - z1p/sin\[Theta]p*vn*vmb);
      Ambmb0S  = Admbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + r1p*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - z1p/sin\[Theta]p*vmb)*vmb;
      FnnIn = RInrp*fnn0;
      FnmbIn = RInrp*fnmb0 - dRInrp*fnmb1;
      FmbmbIn = RInrp*fmbmb0 - dRInrp*fmbmb1 + ddRInrp*fmbmb2;
      dFnnIndr = dRInrp*fnn0 + RInrp*dfnn0dr;
      dFnmbIndr = dRInrp*fnmb0 + RInrp*dfnmb0dr - ddRInrp*fnmb1 - dRInrp*dfnmb1dr;
      dFmbmbIndr = dRInrp*fmbmb0 + RInrp*dfmbmb0dr - ddRInrp*fmbmb1 - dRInrp*dfmbmb1dr + dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr;
      dFnnInd\[Theta] = RInrp*dfnn0d\[Theta];
      dFnmbInd\[Theta] = RInrp*dfnmb0d\[Theta] - dRInrp*dfnmb1d\[Theta];
      dFmbmbInd\[Theta] = RInrp*dfmbmb0d\[Theta] - dRInrp*dfmbmb1d\[Theta] + ddRInrp*dfmbmb2d\[Theta];
      sumPlus  += Total[(\[CapitalSigma]+s*\[CapitalSigma]1)*(s*(
        ((Ann0S + Annt\[Phi]S)*    FnnIn   + AnnrS*  dFnnIndr   + Ann\[Theta]S*  dFnnInd\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)*  FnmbIn  + AnmbrS* dFnmbIndr  + Anmb\[Theta]S* dFnmbInd\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*FmbmbIn + AmbmbrS*dFmbmbIndr + Ambmb\[Theta]S*dFmbmbInd\[Theta])) + 
        (Amnn*FnnIn + Amnmb*FnmbIn + Ammbmb*FmbmbIn))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}];  (* Total over all quadrants *)
      FnnUp = RUprp*fnn0;
      FnmbUp = RUprp*fnmb0 - dRUprp*fnmb1;
      FmbmbUp = RUprp*fmbmb0 - dRUprp*fmbmb1 + ddRUprp*fmbmb2;
      dFnnUpdr = dRUprp*fnn0 + RUprp*dfnn0dr;
      dFnmbUpdr = dRUprp*fnmb0 + RUprp*dfnmb0dr - ddRUprp*fnmb1 - dRUprp*dfnmb1dr;
      dFmbmbUpdr = dRUprp*fmbmb0 + RUprp*dfmbmb0dr - ddRUprp*fmbmb1 - dRUprp*dfmbmb1dr + dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr;
      dFnnUpd\[Theta] = RUprp*dfnn0d\[Theta];
      dFnmbUpd\[Theta] = RUprp*dfnmb0d\[Theta] - dRUprp*dfnmb1d\[Theta];
      dFmbmbUpd\[Theta] = RUprp*dfmbmb0d\[Theta] - dRUprp*dfmbmb1d\[Theta] + ddRUprp*dfmbmb2d\[Theta];
      sumMinus += Total[(\[CapitalSigma]+s*\[CapitalSigma]1)*(s*(
        ((Ann0S + Annt\[Phi]S)*    FnnUp   + AnnrS*  dFnnUpdr   + Ann\[Theta]S*  dFnnUpd\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)*  FnmbUp  + AnmbrS* dFnmbUpdr  + Anmb\[Theta]S* dFnmbUpd\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*FmbmbUp + AmbmbrS*dFmbmbUpdr + Ambmb\[Theta]S*dFmbmbUpd\[Theta])) + 
        (Amnn*FnnUp + Amnmb*FnmbUp + Ammbmb*FmbmbUp))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  CPlus  = 2*Pi*(1-s*\[CapitalGamma]1/\[CapitalGamma])*sumPlus/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Amplitudes *)
  CMinus = 2*Pi*(1-s*\[CapitalGamma]1/\[CapitalGamma])*sumMinus/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus,
      "\[ScriptCapitalH]"->CMinus
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"->(k - \[Omega]*dJzdEn[a,s,{En0+s*En1,Lz0+s*Lz1,Kc0+s*Kc1}] - m*dJzdLz[a,s,{En0+s*En1,Lz0+s*Lz1,Kc0+s*Kc1}])/dJzdK[a,s,{En0+s*En1,Lz0+s*Lz1,Kc0+s*Kc1}]*Abs[CPlus]^2/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->(k - \[Omega]*dJzdEn[a,s,{En0+s*En1,Lz0+s*Lz1,Kc0+s*Kc1}] - m*dJzdLz[a,s,{En0+s*En1,Lz0+s*Lz1,Kc0+s*Kc1}])/dJzdK[a,s,{En0+s*En1,Lz0+s*Lz1,Kc0+s*Kc1}]*\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


(* ::Subsubsection::Closed:: *)
(*Analytical trajectory*)


TeukolskySpinModeAnalytical[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbit_]:=Module[{a,p,e,x,spar,En0,Lz0,K0,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],r,z,Ur,Uz,\[CapitalUpsilon]t,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,
    sumPlus,sumMinus,stepsr,stepsz,zlist,ir,iz,wr,wz,rp,zp,sin\[Theta]p,Urp,Uzp,expr,expz,\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,
    RUprp,dRUprp,ddRUprp,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,
    dfnmb0dr,dfnmb1dr,dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],\[Xi]n,\[Xi]mb,vn,vmb,Sln,Smmb,
    Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,S\[Theta]mb,
    Annt\[Phi]S,AnnrS,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmbt\[Phi]S,Ambmb\[Theta]S,CPlus,CMinus,FnnIn,FnmbIn,FmbmbIn,
    dFnnIndr,dFnmbIndr,dFnmbInd\[Theta],dFmbmbInd\[Theta],FnnUp,FnmbUp,FmbmbUp,
    dFnnUpdr,dFnmbUpdr,dFnmbUpd\[Theta],dFmbmbUpd\[Theta]},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["x"];
  spar = orbit["spar"];
  En0 = orbit["EnTilde"]; (* Geodesic constants of motion *)
  Lz0 = orbit["JzTilde"];
  K0 = orbit["KTilde"];
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = orbit["BLFrequencies"]; (* BL frequencies *)
  r = orbit["Trajectoryg"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z = orbit["Trajectoryg"][[3]];
  \[CapitalUpsilon]t = orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"]; (* Geodesic average rate of change of BL time in Mino time *)
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]z; (* Frequency of mode *)
  If[!(\[Omega]\[Element]Reals), 
  Print["Complex frequency"];
  Return[$Failed]];
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue *)
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega],Method->{"NumericalIntegration",
        "Domain"->{"In"->orbit["rootsTilde"],"Up"->orbit["rootsTilde"]}}];
  \[Lambda] = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalEigenvalue[-2,l,m,a*\[Omega]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);  (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (* Outer horizon *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus  = 0; (* Results of the integration are stored in these variables *)
  sumMinus = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  stepsz = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[stepsz]<>" steps in wz"];
  zlist = {};(* List for functions of \[Theta] *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = r[wr];
    Urp = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+K0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2rp*\[Omega];
    V  = -(K^2+4I*(rp-1)*K)/\[CapitalDelta]+8I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation *)
    RInrp    = R["In"][rp]; (* radial function *)
    dRInrp   = R["In"]'[rp];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    (*dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);*)  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R["Up"][rp];
    dRUprp   = R["Up"]'[rp];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    (*dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);*)
    For[iz = 1, iz <= stepsz/2, iz++,(* integration over wz *)
      wz = N[(iz-1/2)*2Pi/stepsz,Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list in the first step *)
        zp = z[wz];
        Uzp = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(K0-a^2*zp^2)];(* Polar geodesic velocity *)
        expz = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"][wz])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"][wz])+2Pi*k*(iz-1/2)/stepsz)];
        sin\[Theta]p = Sqrt[1-zp^2];
        S = SWSH[ArcCos[zp],0];  (*  Polar function S(\[Theta](z))  *)
        dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]);
        d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta];(* second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[zlist,{zp,Uzp,expz,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}],
        {zp,Uzp,expz,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}=zlist[[iz]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]][[1;;3]];
      {dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a][[2;;-1]];
      \[Xi]n = -\[CapitalDelta]/(2*\[CapitalSigma]);
      \[Xi]mb = I*a*sin\[Theta]p/(Sqrt[2]*\[Zeta]);
      vn = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/(2*\[CapitalSigma])*( 1 + 3*spar*En0/(2*Sqrt[K0]) ) - 3*spar/(2*Sqrt[K0])*\[Xi]n; (* Four-velocity in Kinnersley tetrad *)
      vmb = (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uzp)/(-Sqrt[2]*sin\[Theta]p*\[Zeta])*( 1 + 3*spar*En0/(2*Sqrt[K0]) ) - 3*spar/(2*Sqrt[K0])*\[Xi]mb;
      Sln  = rp/Sqrt[K0]; (* Spin tensor in Kinnersley tetrad *)
      Smmb = I*a*zp/Sqrt[K0];
      rho = 1/\[Zeta]; (* Spin coefficients (with opposite sign than in Teukolsky paper) *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Smmb*(-alpha+Conjugate[beta]);
      Adnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn); (* Dipole terms *)
      Admbmb= (Scd\[Gamma]mbdc*vmb+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Adnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)
        +Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln; (* Functions from Eq. (B3) in [2303.16798] *)
      St\[Phi]mb = (a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      Amnn   = vn^2 + vn*spar/Sqrt[K0]*( \[Xi]n - En0*vn); (* Monopole terms *)
      Amnmb  = vn*vmb + spar/Sqrt[K0]/2*( vn( \[Xi]mb - En0*vmb ) + ( \[Xi]n - En0*vn)*vmb );
      Ammbmb = vmb^2 + vmb*spar/Sqrt[K0]*( \[Xi]mb - En0*vmb );
      Annt\[Phi]S = St\[Phi]n*vn;
      AnnrS  = Srn*vn;
      Anmbt\[Phi]S = (St\[Phi]n*vmb + St\[Phi]mb*vn)/2;
      AnmbrS  = (Srn*vmb)/2;
      Anmb\[Theta]S  = (S\[Theta]mb*vn)/2;
      Ambmbt\[Phi]S = St\[Phi]mb*vmb;
      Ambmb\[Theta]S  = S\[Theta]mb*vmb;
      FnnIn = RInrp*fnn0;
      FnmbIn = RInrp*fnmb0 - dRInrp*fnmb1;
      FmbmbIn = RInrp*fmbmb0 - dRInrp*fmbmb1 + ddRInrp*fmbmb2;
      dFnnIndr = dRInrp*fnn0 + RInrp*dfnn0dr;
      dFnmbIndr = dRInrp*fnmb0 + RInrp*dfnmb0dr - ddRInrp*fnmb1 - dRInrp*dfnmb1dr;
      dFnmbInd\[Theta] = RInrp*dfnmb0d\[Theta] - dRInrp*dfnmb1d\[Theta];
      dFmbmbInd\[Theta] = RInrp*dfmbmb0d\[Theta] - dRInrp*dfmbmb1d\[Theta] + ddRInrp*dfmbmb2d\[Theta];
      sumPlus  += Total[\[CapitalSigma]*( 1 - 3*spar*En0/(2*Sqrt[K0]) )*(spar*(
        ((Adnn + Annt\[Phi]S)*     FnnIn   + AnnrS*  dFnnIndr   +      0*  0       ) +
        ((Adnmb + Anmbt\[Phi]S)*   FnmbIn  + AnmbrS* dFnmbIndr  + Anmb\[Theta]S* dFnmbInd\[Theta]) +
        ((Admbmb + Ambmbt\[Phi]S)* FmbmbIn +      0*0           + Ambmb\[Theta]S*dFmbmbInd\[Theta])) + 
        (Amnn*FnnIn + Amnmb*FnmbIn + Ammbmb*FmbmbIn))*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];  (* Total over all quadrants *)
      FnnUp = RUprp*fnn0;
      FnmbUp = RUprp*fnmb0 - dRUprp*fnmb1;
      FmbmbUp = RUprp*fmbmb0 - dRUprp*fmbmb1 + ddRUprp*fmbmb2;
      dFnnUpdr = dRUprp*fnn0 + RUprp*dfnn0dr;
      dFnmbUpdr = dRUprp*fnmb0 + RUprp*dfnmb0dr - ddRUprp*fnmb1 - dRUprp*dfnmb1dr;
      dFnmbUpd\[Theta] = RUprp*dfnmb0d\[Theta] - dRUprp*dfnmb1d\[Theta];
      dFmbmbUpd\[Theta] = RUprp*dfmbmb0d\[Theta] - dRUprp*dfmbmb1d\[Theta] + ddRUprp*dfmbmb2d\[Theta];
      sumMinus += Total[\[CapitalSigma]*( 1 - 3*spar*En0/(2*Sqrt[K0]) )*(spar*(
        ((Adnn + Annt\[Phi]S)*     FnnUp   + AnnrS*  dFnnUpdr   +      0*  0       ) +
        ((Adnmb + Anmbt\[Phi]S)*   FnmbUp  + AnmbrS* dFnmbUpdr  + Anmb\[Theta]S* dFnmbUpd\[Theta]) +
        ((Admbmb + Ambmbt\[Phi]S)* FmbmbUp +      0*0           + Ambmb\[Theta]S*dFmbmbUpd\[Theta])) + 
        (Amnn*FnnUp + Amnmb*FnmbUp + Ammbmb*FmbmbUp))*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  CPlus  = 2*Pi*sumPlus/(\[CapitalUpsilon]t*W*stepsr*stepsz); (* Amplitudes *)
  CMinus = 2*Pi*sumMinus/(\[CapitalUpsilon]t*W*stepsr*stepsz);
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus,
      "\[ScriptCapitalH]"->CMinus
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "stepsz"->stepsz
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


TeukolskySpinModeAnalyticalNew[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbit_]:=Module[{a,p,e,x,spar,En0,Lz0,K0,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],r,z,\[CapitalUpsilon]t,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,
    sumPlus,sumMinus,stepsr,stepsz,zlist,ir,iz,wr,wz,rp,zp,sin\[Theta]p,Urp,Uzp,expr,expz,\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,
    RUprp,dRUprp,ddRUprp,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,
    dfnmb0dr,dfnmb1dr,dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vn,vmb,
    Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,S\[Theta]mb,
    CPlus,CMinus,FnnIn,FnmbIn,FmbmbIn,
    dFnnIndr,dFnmbIndr,dFnmbInd\[Theta],dFmbmbInd\[Theta],FnnUp,FnmbUp,FmbmbUp,
    dFnnUpdr,dFnmbUpdr,dFnmbUpd\[Theta],dFmbmbUpd\[Theta]},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["x"];
  spar = orbit["spar"];
  En0 = orbit["EnTilde"]; (* Shifted constants of motion *)
  Lz0 = orbit["JzTilde"];
  K0 = orbit["KTilde"];
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = orbit["BLFrequencies"]; (* BL frequencies *)
  r = orbit["Trajectoryg"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z = orbit["Trajectoryg"][[3]];
  \[CapitalUpsilon]t = orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"]; (* Geodesic average rate of change of BL time in Mino time*)
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]z; (* Frequency of mode *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue *)
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega],Method->{"NumericalIntegration",
        "Domain"->{"In"->orbit["rootsTilde"],"Up"->orbit["rootsTilde"]}}];
  \[Lambda] = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalEigenvalue[-2,l,m,a*\[Omega]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);  (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (* Outer horizon *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus  = 0; (* Results of the integration are stored in these variables *)
  sumMinus = 0;
  (* numbers of steps for wr and wz integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  stepsz = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[stepsz]<>" steps in wz"];
  zlist = {};(* List for functions of \[Theta] *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = r[wr];
    Urp = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+K0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2rp*\[Omega];
    V  = -(K^2+4I*(rp-1)*K)/\[CapitalDelta]+8I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation *)
    RInrp    = R["In"][rp]; (* radial function *)
    dRInrp   = R["In"]'[rp];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    (*dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);*)  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R["Up"][rp];
    dRUprp   = R["Up"]'[rp];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    (*dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);*)
    For[iz = 1, iz <= stepsz/2, iz++,(* integration over w\[Theta] *)
      wz = N[(iz-1/2)*2Pi/stepsz,Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list in the first step *)
        zp = z[wz];
        Uzp = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(K0-a^2*zp^2)];(* Polar geodesic velocity *)
        expz = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"][wz])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"][wz])+2Pi*k*(iz-1/2)/stepsz)];
        sin\[Theta]p = Sqrt[1-zp^2];
        S = SWSH[ArcCos[zp],0];  (*  Polar function S(\[Theta](z))  *)
        dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]);
        d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta];(* second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[zlist,{zp,Uzp,expz,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}],
        {zp,Uzp,expz,sin\[Theta]p,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,L2S,dL2Sd\[Theta],L1L2S,dL1L2Sd\[Theta]}=zlist[[iz]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]][[1;;3]];
      {dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a][[2;;-1]];
      vn = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vmb = (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uzp)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
      Adnn  = (rp (-rp \[CapitalDelta]+d\[CapitalDelta] \[CapitalSigma]))/\[CapitalSigma]^2*vn; (* Dipole terms *)
      Admbmb= (Sqrt[2] a zp (a rp (1-zp^2) + I zp \[CapitalSigma]))/(sin\[Theta]p \[Zeta] \[CapitalSigma])*vmb;
      Adnmb = (a (rp (I rp+3 a zp) (1-zp^2) + I zp^2 \[CapitalSigma]))/(2 Sqrt[2] sin\[Theta]p \[Zeta] \[CapitalSigma])*vn + (-a zp (3 I rp+a zp) \[CapitalDelta] + d\[CapitalDelta] rp \[CapitalSigma])/(4*\[CapitalSigma]^2)*vmb;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*rp; (* Functions from Eq. (B3) in [2303.16798] *)
      St\[Phi]mb = (a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*I*a*zp;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*rp;
      S\[Theta]mb = I*a*zp/(Sqrt[2]*\[Zeta]);
      FnnIn = RInrp*fnn0;
      FnmbIn = RInrp*fnmb0 - dRInrp*fnmb1;
      FmbmbIn = RInrp*fmbmb0 - dRInrp*fmbmb1 + ddRInrp*fmbmb2;
      dFnnIndr = dRInrp*fnn0 + RInrp*dfnn0dr;
      dFnmbIndr = dRInrp*fnmb0 + RInrp*dfnmb0dr - ddRInrp*fnmb1 - dRInrp*dfnmb1dr;
      dFnmbInd\[Theta] = RInrp*dfnmb0d\[Theta] - dRInrp*dfnmb1d\[Theta];
      dFmbmbInd\[Theta] = RInrp*dfmbmb0d\[Theta] - dRInrp*dfmbmb1d\[Theta] + ddRInrp*dfmbmb2d\[Theta];
      sumPlus  += Total[\[CapitalSigma]*(spar/Sqrt[K0]*(
        ((Adnn + St\[Phi]n*vn)*     FnnIn   + Srn*vn   *dFnnIndr  +      0*  0       ) +
        ((Adnmb + (St\[Phi]n*vmb + St\[Phi]mb*vn)/2)*   FnmbIn  + Srn*vmb/2*dFnmbIndr + S\[Theta]mb*vn/2*dFnmbInd\[Theta]) +
        ((Admbmb + St\[Phi]mb*vmb)* FmbmbIn +      0*0            + S\[Theta]mb*vmb*dFmbmbInd\[Theta])) + 
        (vn^2*FnnIn + vn*vmb*FnmbIn + vmb^2*FmbmbIn))*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];  (* Total over all quadrants *)
      FnnUp = RUprp*fnn0;
      FnmbUp = RUprp*fnmb0 - dRUprp*fnmb1;
      FmbmbUp = RUprp*fmbmb0 - dRUprp*fmbmb1 + ddRUprp*fmbmb2;
      dFnnUpdr = dRUprp*fnn0 + RUprp*dfnn0dr;
      dFnmbUpdr = dRUprp*fnmb0 + RUprp*dfnmb0dr - ddRUprp*fnmb1 - dRUprp*dfnmb1dr;
      dFnmbUpd\[Theta] = RUprp*dfnmb0d\[Theta] - dRUprp*dfnmb1d\[Theta];
      dFmbmbUpd\[Theta] = RUprp*dfmbmb0d\[Theta] - dRUprp*dfmbmb1d\[Theta] + ddRUprp*dfmbmb2d\[Theta];
      sumMinus += Total[\[CapitalSigma]*(spar/Sqrt[K0]*(
        ((Adnn + St\[Phi]n*vn)*     FnnUp   + Srn*vn*  dFnnUpdr   +      0*  0       ) +
        ((Adnmb + (St\[Phi]n*vmb + St\[Phi]mb*vn)/2)*   FnmbUp  + Srn*vmb/2* dFnmbUpdr  + S\[Theta]mb*vn/2* dFnmbUpd\[Theta]) +
        ((Admbmb + St\[Phi]mb*vmb)* FmbmbUp +      0*0           + S\[Theta]mb*vmb*dFmbmbUpd\[Theta])) + 
        (vn^2*FnnUp + vn*vmb*FnmbUp + vmb^2*FmbmbUp))*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  CPlus  = ( 1 + spar*En0/(2*Sqrt[K0]) )*2*Pi*sumPlus/(\[CapitalUpsilon]t*W*stepsr*stepsz); (* Amplitudes *)
  CMinus = ( 1 + spar*En0/(2*Sqrt[K0]) )*2*Pi*sumMinus/(\[CapitalUpsilon]t*W*stepsr*stepsz);
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus,
      "\[ScriptCapitalH]"->CMinus
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "stepsz"->stepsz
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


TeukolskySpinModeAnalyticalNew2[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbit_]:=Module[{a,p,e,x,spar,
    En0,Lz0,K0,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],r,z,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,sumPlus,sumMinus,stepsr,stepsz,zlist,ir,iz,wr,wz,rp,zp,sin\[Theta]p,K\[Theta],
    Urp,Uzp,expr,expz,\[CapitalDelta],d\[CapitalDelta],K,dK,V,RInrp,dRInrp,ddRInrp,RUprp,dRUprp,ddRUprp,DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,
    \[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,dL2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],vn,vmb,CPlus,CMinus,FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,(*dFnnIndr,dFnmbIndr,dFnmbInd\[Theta],dFmbmbInd\[Theta],*)
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,(*dFnnUpdr,dFnmbUpdr,dFnmbUpd\[Theta],dFmbmbUpd\[Theta],*)(*DFnnIn,*)(*DFnnUp,L0FmbmbUp,*)(*L0FmbmbIn,DFnmbIn,*)(*DFnmbUp,*)(*LFnmbIn,*)(*LFnmbUp,*)
    \[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]z},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["x"];
  spar = orbit["spar"];
  En0 = orbit["EnTilde"]; (* Shifted constants of motion *)
  Lz0 = orbit["JzTilde"];
  K0 = orbit["KTilde"];
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = orbit["BLFrequencies"]; (* BL frequencies *)
  r = orbit["Trajectoryg"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z = orbit["Trajectoryg"][[3]];
  \[CapitalUpsilon]r=orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"];
  \[CapitalUpsilon]z=orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(z\)]\)"];
  \[CapitalUpsilon]\[Phi]=orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"];
  \[CapitalUpsilon]t = orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"]; (* Geodesic average rate of change of BL time in Mino time*)
  \[CapitalUpsilon]tz = orbit["\[CapitalUpsilon]tz"];
  \[CapitalUpsilon]\[Phi]z = orbit["\[CapitalUpsilon]\[Phi]z"];
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]z; (* Frequency of mode *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue *)
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega],Method->{"NumericalIntegration",
        "Domain"->{"In"->orbit["rootsTilde"],"Up"->orbit["rootsTilde"]}}];
  \[Lambda] = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalEigenvalue[-2,l,m,a*\[Omega]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);  (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (* Outer horizon *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus  = 0; (* Results of the integration are stored in these variables *)
  sumMinus = 0;
  (* numbers of steps for wr and wz integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  stepsz = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[stepsz]<>" steps in wz"];
  zlist = {};(* List for functions of \[Theta] *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = r[wr];
    Urp = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+K0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2rp*\[Omega];
    V  = -(K^2+4I*(rp-1)*K)/\[CapitalDelta]+8I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
    RInrp    = R["In"][rp]; (* radial function *)
    dRInrp   = R["In"]'[rp];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    RUprp    = R["Up"][rp];
    dRUprp   = R["Up"]'[rp];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    DRIn = dRInrp - I*K/\[CapitalDelta]*RInrp;
    dDRIndr = ddRInrp - I*K/\[CapitalDelta]*dRInrp - I*( dK*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RInrp;
    DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
    DRUp = dRUprp - I*K/\[CapitalDelta]*RUprp;
    dDRUpdr = ddRUprp - I*K/\[CapitalDelta]*dRUprp - I*( dK*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUprp;
    DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
    For[iz = 1, iz <= stepsz/2, iz++,(* integration over w\[Theta] *)
      wz = N[(iz-1/2)*2Pi/stepsz,Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list in the first step *)
        zp = z[wz];
        Uzp = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(K0-a^2*zp^2)];(* Polar geodesic velocity *)
        expz = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tz"][wz])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]z"][wz])+2Pi*k*(iz-1/2)/stepsz)];
        sin\[Theta]p = Sqrt[1-zp^2];
        K\[Theta] = a*\[Omega]*(1-zp^2) - m;
        S = SWSH[ArcCos[zp],0];  (*  Polar function S(\[Theta](z))  *)
        dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]);
        d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta];(* second derivative of S from Teukolsky equation *)
        L2S = dSd\[Theta] + (K\[Theta] + 2 zp)*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        dL2Sd\[Theta] = d2Sd\[Theta]2 + (K\[Theta] + 2 zp)*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*S;
        L1L2S = dL2Sd\[Theta] + (K\[Theta] + zp)*L2S/sin\[Theta]p; 
        AppendTo[zlist,{zp,Uzp,expz,sin\[Theta]p,K\[Theta],S,dSd\[Theta],d2Sd\[Theta]2,L2S,dL2Sd\[Theta],L1L2S}],
        {zp,Uzp,expz,sin\[Theta]p,K\[Theta],S,dSd\[Theta],d2Sd\[Theta]2,L2S,dL2Sd\[Theta],L1L2S}=zlist[[iz]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      vn = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/\[CapitalDelta]; (* Four-velocity in Kinnersley tetrad *)
      vmb = (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uzp)/sin\[Theta]p^2;
      {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RInrp,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
      sumPlus  += Total[(
          vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn + spar/Sqrt[K0]*( vn*GnIn + vmb*GmbIn )
         )*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];  (* Total over all quadrants *)
      {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RUprp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
      sumMinus += Total[(
          vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp + spar/Sqrt[K0]*( vn*GnUp + vmb*GmbUp )
         )*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  CPlus  = ( 1 + spar*En0/(2*Sqrt[K0]) )*2*Pi*sumPlus/(\[CapitalUpsilon]t*W*stepsr*stepsz); (* Amplitudes *)
  CMinus = ( 1 + spar*En0/(2*Sqrt[K0]) )*2*Pi*sumMinus/(\[CapitalUpsilon]t*W*stepsr*stepsz);
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus,
      "\[ScriptCapitalH]"->CMinus
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"->(1 + spar*En0/(2*Sqrt[K0]))*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*Abs[CPlus]^2/(2Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->(1 + spar*En0/(2*Sqrt[K0]))*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*\[Alpha]*Abs[CMinus]^2/(2Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "stepsz"->stepsz
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


(* ::Subsubsection::Closed:: *)
(*Corrections*)


TeukolskySpinModeCorrectionNum[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbitCorrection_,\[Delta]\[Omega]_]:=Module[{h1,h2,h3,a,p,e,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,correction,r,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,
    SWSH,SWSHplus,SWSHminus,R,Rplus,Rminus,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    correctionp,lists,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,
    RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,\[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_]:=1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  h3[r_,z_]:=-((2 a r z)/(Sqrt[Kc] (r^2+a^2 z^2)^2));
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  \[ScriptCapitalI] = orbitCorrection["\[ScriptCapitalI]"];
  En = orbitCorrection["Ehat"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzhat"];
  Kc = orbitCorrection["Khat"];
  En1 = orbitCorrection["ES"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["LS"];
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Geodesic BL frequencies *)
  {\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  r = orbitCorrection["TrajectoryGeo"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma]  = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega]  = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + n*\[CapitalOmega]r1 + k*\[CapitalOmega]\[Theta]1;
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue for geodesic frequency and for numerical derivative *)
  SWSHplus = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*(\[Omega]+\[Delta]\[Omega])];
  SWSHminus = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*(\[Omega]-\[Delta]\[Omega])];
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega],Method->{"NumericalIntegration","Domain"->{"In"->{p/(1+e),p/(1-e)},"Up"->{p/(1+e),p/(1-e)}}}];
  Rplus = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega]+\[Delta]\[Omega],Method->{"NumericalIntegration","Domain"->{"In"->{p/(1+e),p/(1-e)},"Up"->{p/(1+e),p/(1-e)}}}];
  Rminus = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega]-\[Delta]\[Omega],Method->{"NumericalIntegration","Domain"->{"In"->{p/(1+e),p/(1-e)},"Up"->{p/(1+e),p/(1-e)}}}];
  \[Lambda] = R["In"]["Eigenvalue"];
  \[Lambda]1 = (Rplus["In"]["Eigenvalue"] - Rminus["In"]["Eigenvalue"])/(2*\[Delta]\[Omega])*\[Omega]1; (* Linear correction to the eigenvalue *)
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  sumPlus0 = sumPlus1 = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in w\[Theta]"];
  \[Theta]list = {};(* List for functions of \[Theta] *)
  correctionp = Table[correction[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]],N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]]],{ir,1,stepsr/2},{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points in one quadrant in w_r and w_\[Theta] *)
  For[ir = 1, ir <= stepsr/2, ir++,(* integration over w_r *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]];
    rp = r[wr];
    Ur = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];(* Geodesic radial velocity *)
    expr = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    (*expr1 = {1,-1,1,-1}*I*\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr];(* Linear correction to the exponential term due to change of frequencies *)*)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    K1  = (rp^2+a^2)*\[Omega]1;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2*rp*\[Omega];
    dK1 = 2*rp*\[Omega]1;
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    V1  = -(2*K*K1 + 4I*(rp-1)*K1)/\[CapitalDelta] + 8*I*\[Omega]1*rp + \[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RInrp    = R["In"][rp]; (* radial function *)
    dRInrp   = R["In"]'[rp];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    RInrp1   = (Rplus["In"][rp]-Rminus["In"][rp])/(2*\[Delta]\[Omega])*\[Omega]1; (* Linear part of the radial function *)
    dRInrp1  = (Rplus["In"]'[rp]-Rminus["In"]'[rp])/(2*\[Delta]\[Omega])*\[Omega]1;
    ddRInrp1 = (V*RInrp1+V1*RInrp+d\[CapitalDelta]*dRInrp1)/\[CapitalDelta];
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R["Up"][rp];
    dRUprp   = R["Up"]'[rp];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    RUprp1   = (Rplus["Up"][rp]-Rminus["Up"][rp])/(2*\[Delta]\[Omega])*\[Omega]1;
    dRUprp1  = (Rplus["Up"]'[rp]-Rminus["Up"]'[rp])/(2*\[Delta]\[Omega])*\[Omega]1;
    ddRUprp1 = (V*RUprp1+V1*RUprp+d\[CapitalDelta]*dRUprp1)/\[CapitalDelta];
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++, (* integration over w_\[Theta] *)
      w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta], Precision[{a,p,e,\[ScriptCapitalI]}]];
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        zp = z[w\[Theta]];
        Uz = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2];
        S = SWSH[ArcCos[zp],0];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
        S1 = (SWSHplus[ArcCos[zp],0]-SWSHminus[ArcCos[zp],0])/(2*\[Delta]\[Omega])*\[Omega]1;  (*  Linear part of S(\[Theta](z))  *)
        dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]); (* First derivative of S wrt \[Theta] *)
        dS1d\[Theta] = ((D[SWSHplus[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp])-(D[SWSHminus[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]))/(2*\[Delta]\[Omega])*\[Omega]1; (* Linear part of first derivative of S wrt \[Theta] *)
        d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta]; (* second derivative of S from Teukolsky equation *)
        d2S1d\[Theta]2 = (-(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S1-(-a^2*2*\[Omega]*\[Omega]1*(1-zp^2)+4*a*zp*\[Omega]1+2*m*a*\[Omega]1+\[Lambda]1)*S-zp/sin\[Theta]p*dS1d\[Theta]);(* Linear part of second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        L2S1 = dS1d\[Theta]-(S1 (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2) \[Omega]1))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2) \[Omega]1))/sin\[Theta]p + 
                 S1 (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) \[Omega]1 - a^2 (-1+zp^2) 2*\[Omega]*\[Omega]1);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}];,
        {zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,K1,dK1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
      vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
      Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
      Snmb = Conjugate[Snm];
      Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
      Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
      Amnn   = vn^2;
      Amnmb  = vn*vmb;
      Ammbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Adnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Admbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Adnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      rp1 =  {correctionp[[ir,i\[Theta]]]["rS"],  correctionp[[ir,-i\[Theta]]]["rS"],  correctionp[[ir,-i\[Theta]]]["rS"],  correctionp[[ir,i\[Theta]]]["rS"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
      zp1 =  {correctionp[[ir,i\[Theta]]]["zS"], -correctionp[[ir,-i\[Theta]]]["zS"], -correctionp[[ir,-i\[Theta]]]["zS"],  correctionp[[ir,i\[Theta]]]["zS"]};
      Urp1 = {correctionp[[ir,i\[Theta]]]["UrS"],-correctionp[[ir,-i\[Theta]]]["UrS"], correctionp[[ir,-i\[Theta]]]["UrS"],-correctionp[[ir,i\[Theta]]]["UrS"]};
      Uzp1 = {correctionp[[ir,i\[Theta]]]["UzS"], correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,i\[Theta]]]["UzS"]};
      \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
      exp1 = I*{(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
               -(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
                (\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
               -(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      vn1  = -( 2*rp*rp1*En - ((rp^2+a^2)*En - a*Lz + Ur)/\[CapitalSigma]*\[CapitalSigma]1 + 
               (rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz) + Urp1)/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = ( I*2*a*zp*zp1*En - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
               -I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz)) + Uzp1)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Adnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Adnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Admbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RInrp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRInrp + Ammbmb*fmbmb2*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RUprp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRUprp + Ammbmb*fmbmb2*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*RInrp1 - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*dRInrp1 + Ammbmb*fmbmb2*ddRInrp1
          + (Amnn*fnn01 + Amnmb*fnmb01 + Ammbmb*fmbmb01)*RInrp - (Amnmb*fnmb11 + Ammbmb*fmbmb11)*dRInrp + Ammbmb*fmbmb21*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*RUprp1 - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*dRUprp1 + Ammbmb*fmbmb2*ddRUprp1
          + (Amnn*fnn01 + Amnmb*fnmb01 + Ammbmb*fmbmb01)*RUprp - (Amnmb*fnmb11 + Ammbmb*fmbmb11)*dRUprp + Ammbmb*fmbmb21*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/\[CapitalDelta]; (* Invariant Wronskian *)
  W1 = (RInrp1*dRUprp + RInrp*dRUprp1 - dRInrp1*RUprp - dRInrp*RUprp1)/\[CapitalDelta]; (* Linear part of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - W1/W*sumPlus0 )/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - W1/W*sumMinus0)/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  Association[
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  ] (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


Options[TeukolskySpinModeCorrection] = {WorkingPrecision->32};
TeukolskySpinModeCorrection[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbitCorrection_,{angparNew_,RCorrection_},OptionsPattern[]]:=Module[{h1,h2,h3,a,p,e,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,correction,
    r,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,SWSH,dSWSHd\[Omega],R,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    correctionp,lists,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,
    RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,\[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1,\[Psi],d\[Lambda]d\[Omega]},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_]:=1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+a^7 En^2 z^4 (-1+z^2)
               +a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))
               +a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  h3[r_,z_]:=-((2 a r z)/(Sqrt[Kc] (r^2+a^2 z^2)^2));
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  \[ScriptCapitalI] = orbitCorrection["\[ScriptCapitalI]"];
  En = orbitCorrection["Ehat"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzhat"];
  Kc = orbitCorrection["Khat"];
  En1 = orbitCorrection["ES"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["LS"];
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Geodesic BL frequencies *)
  {\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  r = orbitCorrection["TrajectoryGeo"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma]  = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega]  = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + n*\[CapitalOmega]r1 + k*\[CapitalOmega]\[Theta]1;
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,
                                     SetPrecision[a, OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega], OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                         "precODE" -> OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  R = RCorrection[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],
                         SetPrecision[\[Omega],OptionValue[WorkingPrecision]+5],1,
                         SetPrecision[\[Lambda],OptionValue[WorkingPrecision]+5],
                         SetPrecision[d\[Lambda]d\[Omega],OptionValue[WorkingPrecision]+5],e,p,"precODE"->OptionValue[WorkingPrecision]];
  \[Lambda]1 = d\[Lambda]d\[Omega]*\[Omega]1;
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  sumPlus0 = sumPlus1 = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and w\[Theta] integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in w\[Theta]"];
  \[Theta]list = {};(* List for functions of \[Theta] *)
  correctionp = Table[correction[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]],N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]]],{ir,1,stepsr/2},{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points in one quadrant in w_r and w_\[Theta] *)
  For[ir = 1, ir <= stepsr/2, ir++,(* integration over w_r *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]];
    rp = r[wr];
    Ur = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En-a*Lz)^2-(rp^2-2rp+a^2)*(rp^2+Kc)];(* Geodesic radial velocity *)
    expr = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Geodesic exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    K1  = (rp^2+a^2)*\[Omega]1;
    d\[CapitalDelta] = 2*(rp-1);
    dK = 2*rp*\[Omega];
    dK1 = 2*rp*\[Omega]1;
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    V1  = -(2*K*K1 + 4I*(rp-1)*K1)/\[CapitalDelta] + 8*I*\[Omega]1*rp + \[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
    dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
    RInrp    = R[rp]["R0"]["In"]["R"]; (* radial function *)
    dRInrp   = R[rp]["R0"]["In"]["dR"];
    ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    RInrp1   = R[rp]["R1"]["In"]["R"]*\[Omega]1; (* Linear part of the radial function *)
    dRInrp1  = R[rp]["R1"]["In"]["dR"]*\[Omega]1;
    ddRInrp1 = (V*RInrp1+V1*RInrp+d\[CapitalDelta]*dRInrp1)/\[CapitalDelta];
    dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
    RUprp    = R[rp]["R0"]["Up"]["R"];
    dRUprp   = R[rp]["R0"]["Up"]["dR"];
    ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
    RUprp1   = R[rp]["R1"]["Up"]["R"]*\[Omega]1;
    dRUprp1  = R[rp]["R1"]["Up"]["dR"]*\[Omega]1;
    ddRUprp1 = (V*RUprp1+V1*RUprp+d\[CapitalDelta]*dRUprp1)/\[CapitalDelta];
    dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
    For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++, (* integration over w_\[Theta] *)
      w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta], Precision[{a,p,e,\[ScriptCapitalI]}]];
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        zp = z[w\[Theta]];
        Uz = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
        exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
        sin\[Theta]p = Sqrt[1-zp^2];
        S = SWSH[ArcCos[zp]][[1]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
        S1 = dSWSHd\[Omega][ArcCos[zp]][[1]]*\[Omega]1;  (*  Linear part of S(\[Theta](z))  *)
        dSd\[Theta] = SWSH[ArcCos[zp]][[2]]; (* First derivative of S wrt \[Theta] *)
        dS1d\[Theta] = dSWSHd\[Omega][ArcCos[zp]][[2]]*\[Omega]1; (* Linear part of first derivative of S wrt \[Theta] *)
        d2Sd\[Theta]2 = -(zp/sin\[Theta]p)*dSd\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S; (* second derivative of S from Teukolsky equation *)
        d2S1d\[Theta]2 = (-(zp/sin\[Theta]p)*dS1d\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S1 - (-2*a^2*\[Omega]*(1-zp^2) + 4*a*zp + 2*m*a + d\[Lambda]d\[Omega])*\[Omega]1*S);(* Linear part of second derivative of S from Teukolsky equation *)
        d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
                 -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
        L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        L2S1 = dS1d\[Theta]-(S1 (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2) \[Omega]1))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
        dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
        L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
        L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2) \[Omega]1))/sin\[Theta]p + 
                 S1 (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) \[Omega]1 - a^2 (-1+zp^2) 2*\[Omega]*\[Omega]1);  
        dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
                   1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
        AppendTo[\[Theta]list,{zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}];,
        {zp,Uz,exp\[Theta],exp\[Theta]1,sin\[Theta]p,S,S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,L2S,L2S1,dL2Sd\[Theta],L1L2S,L1L2S1,dL1L2Sd\[Theta]}=\[Theta]list[[i\[Theta]]];
      ];
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      \[CapitalSigma] = rp^2+a^2*zp^2;
      {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
      {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
      {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
      {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,K1,dK1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
      vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
      vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
      vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
      vmb = Conjugate[vm];
      Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
      Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
      Snmb = Conjugate[Snm];
      Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
      Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
      Amnn   = vn^2;
      Amnmb  = vn*vmb;
      Ammbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Adnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Admbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Adnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      rp1 =  {correctionp[[ir,i\[Theta]]]["rS"],  correctionp[[ir,-i\[Theta]]]["rS"],  correctionp[[ir,-i\[Theta]]]["rS"],  correctionp[[ir,i\[Theta]]]["rS"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
      zp1 =  {correctionp[[ir,i\[Theta]]]["zS"], -correctionp[[ir,-i\[Theta]]]["zS"], -correctionp[[ir,-i\[Theta]]]["zS"],  correctionp[[ir,i\[Theta]]]["zS"]};
      Urp1 = {correctionp[[ir,i\[Theta]]]["UrS"],-correctionp[[ir,-i\[Theta]]]["UrS"], correctionp[[ir,-i\[Theta]]]["UrS"],-correctionp[[ir,i\[Theta]]]["UrS"]};
      Uzp1 = {correctionp[[ir,i\[Theta]]]["UzS"], correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,-i\[Theta]]]["UzS"],-correctionp[[ir,i\[Theta]]]["UzS"]};
      \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
      exp1 = I*{(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
               -(\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
                (\[Omega]*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir,-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
               -(\[Omega]*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ir, i\[Theta]]]["\[CapitalDelta]\[Phi]S"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr]-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
      vn1  = -( 2*rp*rp1*En - ((rp^2+a^2)*En - a*Lz + Ur)/\[CapitalSigma]*\[CapitalSigma]1 + 
               (rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz) + Urp1)/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
      vmb1 = ( I*2*a*zp*zp1*En - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
               -I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp]+h3[rp,zp]*Ur*Uz)) + Uzp1)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Adnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Adnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Admbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RInrp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRInrp + Ammbmb*fmbmb2*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RUprp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRUprp + Ammbmb*fmbmb2*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*RInrp1 - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*dRInrp1 + Ammbmb*fmbmb2*ddRInrp1
          + (Amnn*fnn01 + Amnmb*fnmb01 + Ammbmb*fmbmb01)*RInrp - (Amnmb*fnmb11 + Ammbmb*fmbmb11)*dRInrp + Ammbmb*fmbmb21*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*RUprp1 - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*dRUprp1 + Ammbmb*fmbmb2*ddRUprp1
          + (Amnn*fnn01 + Amnmb*fnmb01 + Ammbmb*fmbmb01)*RUprp - (Amnmb*fnmb11 + Ammbmb*fmbmb11)*dRUprp + Ammbmb*fmbmb21*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
    ];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/\[CapitalDelta]; (* Invariant Wronskian *)
  W1 = (RInrp1*dRUprp + RInrp*dRUprp1 - dRInrp1*RUprp - dRInrp*RUprp1)/\[CapitalDelta]; (* Linear part of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - W1/W*sumPlus0 )/(\[CapitalGamma]*W*stepsr*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - W1/W*sumMinus0)/(\[CapitalGamma]*W*stepsr*steps\[Theta]);
  Association[
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->steps\[Theta]
  ] (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


Options[TeukolskySpinModeCorrectionAnalytical] = {WorkingPrecision->32};
TeukolskySpinModeCorrectionAnalytical[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbit_,{angparNew_,RCorrection_},OptionsPattern[]]:=Module[{
    a,p,e,x,En0,Lz0,K0,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],r,z,\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t1,\[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]z,\[CapitalUpsilon]\[Tau]z,\[Omega],\[Omega]1,
    SWSH,dSWSHd\[Omega],R,\[Lambda],d\[Lambda]d\[Omega],\[ScriptCapitalC]2,d\[ScriptCapitalC]2d\[Omega],rplus,P,\[Epsilon],\[Alpha],d\[Alpha]d\[Omega],
    sumPlus0,sumPlus1,sumMinus0,sumMinus1,stepsr,stepsz,zlist,ir,iz,wr,wz,
    rp,Urp,expr,expr1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],V,dVd\[Omega],
    Rrp,RIn,dRIndr,d2RIndr2,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,dDRInd\[Omega],d2DRIndrd\[Omega],dDDRInd\[Omega],dDRUpd\[Omega],d2DRUpdrd\[Omega],dDDRUpd\[Omega],
    zp,\[Theta]p,Uzp,expz,expz1,sin\[Theta]p,K\[Theta],dK\[Theta]d\[Omega],S,dSd\[Theta],d2Sd\[Theta]2,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],
    L2S,dL2Sd\[Theta],L1L2S,dL2Sd\[Omega],d2L2Sd\[Theta]d\[Omega],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,vn,vmb,
    FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega],
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega],
    W,W1,CPlus0,CPlus1,CMinus0,CMinus1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["Inclination"];
  En0 = orbit["Energy"]; (* Shifted constants of motion *)
  Lz0 = orbit["AngularMomentum"];
  K0 = orbit["CarterConstant"] + (Lz0 - a*En0)^2;
  {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[orbit["Frequencies"]];
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t; (* BL frequencies *)
  r = orbit["Trajectory"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[qz_] := Cos[orbit["Trajectory"][[3]][qz]];
  \[CapitalUpsilon]\[Tau] = orbit["ProperTimeFrequency"];
  \[CapitalUpsilon]t1 = -3*\[CapitalUpsilon]\[Tau]/(2*Sqrt[K0]);
  \[CapitalUpsilon]tz = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencyt\[Theta][a,p,e,x,{En0,Lz0,K0 - (Lz0 - a*En0)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,p,e,x]] + a*Lz0;
  \[CapitalUpsilon]\[Phi]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Phi]\[Theta][a,p,e,x,{En0,Lz0,K0 - (Lz0 - a*En0)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,p,e,x]] - a*En0;
  \[CapitalUpsilon]\[Tau]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Tau]z[a,p,e,x];
  Print["\[CapitalUpsilon]\[Tau]z = "<>ToString[\[CapitalUpsilon]\[Tau]z]];
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]z; (* Frequency of mode *)
  \[Omega]1 = 3*\[CapitalUpsilon]\[Tau]*\[Omega]/(2*Sqrt[K0]*\[CapitalUpsilon]t);
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,
                                     SetPrecision[a, OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega], OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                         "precODE" -> OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  R = RCorrection[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],
                         SetPrecision[\[Omega],OptionValue[WorkingPrecision]+5],1,
                         SetPrecision[\[Lambda],OptionValue[WorkingPrecision]+5],
                         SetPrecision[d\[Lambda]d\[Omega],OptionValue[WorkingPrecision]+5],e,p,"precODE"->OptionValue[WorkingPrecision]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  d\[ScriptCapitalC]2d\[Omega] = 4 \[Lambda]^3 d\[Lambda]d\[Omega]+4 \[Lambda]^2 (3 d\[Lambda]d\[Omega] + 10 a (m-2 a \[Omega]))+8 \[Lambda] (d\[Lambda]d\[Omega] (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega])) + 
        48 \[Omega] (a m d\[Lambda]d\[Omega]+6-18 a^3 m \[Omega]+12 a^4 \[Omega]^2+a^2 (d\[Lambda]d\[Omega] \[Omega]+6 m^2));  (* \[Omega] derivative of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  d\[Alpha]d\[Omega] = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*d\[ScriptCapitalC]2d\[Omega] + 256*(2*rplus)^5*((P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2)/\[ScriptCapitalC]2; (* \[Omega] derivative of the constant for horizon fluxes *)
  sumPlus0  = sumPlus1  = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and wz integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  stepsz = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[stepsz]<>" steps in wz"];
  zlist = {};(* List for functions of \[Theta] *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = r[wr];
    Urp = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+K0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    expr1 = {1,-1,1,-1}*3*I*\[Omega]/(2*Sqrt[K0])*(\[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr] - orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]r"][wr]);(* Linear part of the exponential term *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    dKdr = 2*rp*\[Omega];
    dKd\[Omega]  = (rp^2+a^2);
    d2Kdrd\[Omega] = 2*rp;
    d\[CapitalDelta] = 2*(rp-1);
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
    Rrp = R[rp];
    (* In solutions of the radial equation and their r and \[Omega] derivatives *)
    RIn    = Rrp["R0"]["In"]["R"]; 
    dRIndr   = Rrp["R0"]["In"]["dR"];
    d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dRInd\[Omega]   = Rrp["R1"]["In"]["R"]; (* \[Omega]-derivatives of the radial function *)
    d2RIndrd\[Omega]  = Rrp["R1"]["In"]["dR"];
    d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
    (* Operators D and their r and \[Omega] derivatives *)
    DRIn = dRIndr - I*K/\[CapitalDelta]*RIn;
    dDRInd\[Omega] = d2RIndrd\[Omega] - I*K/\[CapitalDelta]*dRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RIn;
    dDRIndr = d2RIndr2 - I*K/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
    d2DRIndrd\[Omega] = d3RIndr2d\[Omega] - I*K/\[CapitalDelta]*d2RIndrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRInd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
    DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
    dDDRInd\[Omega] = d2DRIndrd\[Omega] - I*K/\[CapitalDelta]*dDRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRIn;
    (* Up solutions of the radial equation and their r and \[Omega] derivatives *)
    RUp    = Rrp["R0"]["Up"]["R"]; 
    dRUpdr   = Rrp["R0"]["Up"]["dR"];
    d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dRUpd\[Omega]   = Rrp["R1"]["Up"]["R"]; (* \[Omega]-derivative of the radial function *)
    d2RUpdrd\[Omega]  = Rrp["R1"]["Up"]["dR"];
    d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
    (* Operators D and their r and \[Omega] derivatives *)
    DRUp = dRUpdr - I*K/\[CapitalDelta]*RUp;
    dDRUpd\[Omega] = d2RUpdrd\[Omega] - I*K/\[CapitalDelta]*dRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RUp;
    dDRUpdr = d2RUpdr2 - I*K/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
    d2DRUpdrd\[Omega] = d3RUpdr2d\[Omega] - I*K/\[CapitalDelta]*d2RUpdrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRUpd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
    DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
    dDDRUpd\[Omega] = d2DRUpdrd\[Omega] - I*K/\[CapitalDelta]*dDRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRUp;
    (* integration over w\[Theta] *)
    For[iz = 1, iz <= stepsz/2, iz++,
      wz = N[(iz-1/2)*2Pi/stepsz,Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list in the first step *)
        zp = z[wz];
        \[Theta]p = ArcCos[zp];
        Uzp = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(K0-a^2*zp^2)];(* Polar geodesic velocity *)
        expz = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][wz])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][wz])+2Pi*k*(iz-1/2)/stepsz)];(* Exponential term with geodesic \[CapitalDelta]tz and \[CapitalDelta]\[Phi]z *)
        expz1 = {1,1,-1,-1}*3*I*\[Omega]/(2*Sqrt[K0])*(\[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][wz] - orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]\[Theta]"][wz]);(* Linear part of the exponential term *)
        sin\[Theta]p = Sqrt[1-zp^2];
        K\[Theta] = a*\[Omega]*(1-zp^2) - m;
        dK\[Theta]d\[Omega] = a*(1-zp^2);
        (* Spin-weighted spheroidal harmonics and their \[Theta]-derivatives *)
        {S,dSd\[Theta]}=SWSH[\[Theta]p]; 
        d2Sd\[Theta]2 = -(zp/sin\[Theta]p)*dSd\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
        (* \[Omega]-derivatives of S *)
        {dSd\[Omega],d2Sd\[Theta]d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
        d3Sd\[Theta]2d\[Omega] = -(zp/sin\[Theta]p)*d2Sd\[Theta]d\[Omega] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega]*(1-zp^2) + 4*a*zp + 2*m*a + d\[Lambda]d\[Omega])*S;
        (* Operators L and their \[Omega]-derivatives *)
        L2S = dSd\[Theta] + (K\[Theta] + 2 zp)*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        dL2Sd\[Theta] = d2Sd\[Theta]2 + (K\[Theta] + 2 zp)*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*S;
        L1L2S = dL2Sd\[Theta] + (K\[Theta] + zp)*L2S/sin\[Theta]p; 
        dL2Sd\[Omega] = d2Sd\[Theta]d\[Omega] + (K\[Theta] + 2 zp)*dSd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        d2L2Sd\[Theta]d\[Omega] = d3Sd\[Theta]2d\[Omega] + (K\[Theta] + 2 zp)*d2Sd\[Theta]d\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*dSd\[Omega] + (2*a*zp - dK\[Theta]d\[Omega]/sin\[Theta]p^2*zp)*S;
        dL1L2Sd\[Omega] = d2L2Sd\[Theta]d\[Omega] + (K\[Theta] + zp)*dL2Sd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*L2S/sin\[Theta]p; 
        AppendTo[zlist,{zp,Uzp,expz,expz1,sin\[Theta]p,S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]}],
        {zp,Uzp,expz,expz1,sin\[Theta]p,S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]}=zlist[[iz]];
      ];
      (* Integration f functions of r and \[Theta] *)
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      vn = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
      vmb = (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uzp)/sin\[Theta]p^2;
      (* In solution *)
      {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RIn,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
      {dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega]}   = dFabd\[Omega][\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RIn,DRIn,DDRIn,dRInd\[Omega],dDRInd\[Omega],dDDRInd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
      (* Summing the geodesic ampltude over all quadrants *)
      sumPlus0  += Total[(vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];
      (* Summing the linear part of the ampltude over all quadrants *)
      sumPlus1  += Total[(
          (vn*( GnIn ) + vmb*( GmbIn ))/Sqrt[K0] +
          (vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*(expr1 + expz1) +
          (vn*vn*dFnnInd\[Omega] + vn*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
         )*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}]; 
      (* Up solution *)
      {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RUp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
      {dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega]}   = dFabd\[Omega][\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RUp,DRUp,DDRUp,dRUpd\[Omega],dDRUpd\[Omega],dDDRUpd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
      (* Summing the geodesic ampltude over all quadrants *)
      sumMinus0  += Total[(vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];
      (* Summing the linear part of the ampltude over all quadrants *)
      sumMinus1  += Total[(
          (vn*( GnUp ) + vmb*( GmbUp ))/Sqrt[K0] +
          (vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*(expr1 + expz1) +
          (vn*vn*dFnnUpd\[Omega] + vn*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
         )*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}]; 
    ];
  ];
  W = (RIn*dRUpdr - dRIndr*RUp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  W1 = (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/(rp^2-2*rp+a^2)*\[Omega]1; (* Invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalUpsilon]t*W*stepsr*stepsz); (* Amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalUpsilon]t*W*stepsr*stepsz);
  CPlus1  = 2*Pi*(sumPlus1)/(\[CapitalUpsilon]t*W*stepsr*stepsz) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CPlus0; (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1)/(\[CapitalUpsilon]t*W*stepsr*stepsz) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CMinus0;
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"->(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*Abs[CPlus0]^2/(2Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*\[Alpha]*Abs[CMinus0]^2/(2Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"-> (En0*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z) - a*( m - a*\[Omega] ) - 3*\[Omega]*( -\[CapitalUpsilon]\[Tau]z + \[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*\[CapitalUpsilon]tz ))/(2*Sqrt[K0])*Abs[CPlus0]^2/(2Pi*\[Omega]^3)+(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(2Pi*\[Omega]^3),
        "\[ScriptCapitalH]"-> (En0*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z) - a*( m - a*\[Omega] ) - 3*\[Omega]*( -\[CapitalUpsilon]\[Tau]z + \[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*\[CapitalUpsilon]tz ))/(2*Sqrt[K0])*\[Alpha]*Abs[CMinus0]^2/(2Pi*\[Omega]^3)+(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(2Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->stepsz
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


Options[TeukolskySpinModeCorrectionAnalyticalNew] = {WorkingPrecision->32};
TeukolskySpinModeCorrectionAnalyticalNew[l_?IntegerQ,m_?IntegerQ,n_?IntegerQ,k_?IntegerQ,orbit_,{angparNew_,RCorrection_},OptionsPattern[]]:=Module[{
    a,p,e,x,En0,Lz0,K0,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],r,z,\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t1,\[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]z,\[CapitalUpsilon]\[Tau]z,\[Omega],\[Omega]1,
    SWSH,dSWSHd\[Omega],R,\[Lambda],d\[Lambda]d\[Omega],\[ScriptCapitalC]2,d\[ScriptCapitalC]2d\[Omega],rplus,P,\[Epsilon],\[Alpha],d\[Alpha]d\[Omega],
    sumPlus0,sumPlus1,sumMinus0,sumMinus1,stepsr,stepsz,zlist,ir,iz,wr,wz,
    rp,Urp,expr,expr1,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],V,dVd\[Omega],
    RIn,dRIndr,d2RIndr2,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,dDRInd\[Omega],d2DRIndrd\[Omega],dDDRInd\[Omega],dDRUpd\[Omega],d2DRUpdrd\[Omega],dDDRUpd\[Omega],
    zp,\[Theta]p,Uzp,expz,expz1,sin\[Theta]p,K\[Theta],dK\[Theta]d\[Omega],S,dSd\[Theta],d2Sd\[Theta]2,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],
    L2S,dL2Sd\[Theta],L1L2S,dL2Sd\[Omega],d2L2Sd\[Theta]d\[Omega],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,vn,vmb,
    FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega],
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega],
    W,W1,CPlus0,CPlus1,CMinus0,CMinus1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["Inclination"];
  En0 = orbit["Energy"]; (* Shifted constants of motion *)
  Lz0 = orbit["AngularMomentum"];
  K0 = orbit["CarterConstant"] + (Lz0 - a*En0)^2;
  {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[orbit["Frequencies"]];
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t; (* BL frequencies *)
  r = orbit["Trajectory"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[qz_] := Cos[orbit["Trajectory"][[3]][qz]];
  \[CapitalUpsilon]\[Tau] = orbit["ProperTimeFrequency"];
  \[CapitalUpsilon]t1 = -3*\[CapitalUpsilon]\[Tau]/(2*Sqrt[K0]);
  \[CapitalUpsilon]tz = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencyt\[Theta][a,p,e,x,{En0,Lz0,K0 - (Lz0 - a*En0)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,p,e,x]] + a*Lz0;
  \[CapitalUpsilon]\[Phi]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Phi]\[Theta][a,p,e,x,{En0,Lz0,K0 - (Lz0 - a*En0)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,p,e,x]] - a*En0;
  \[CapitalUpsilon]\[Tau]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Tau]z[a,p,e,x];
  Print["\[CapitalUpsilon]\[Tau]z = "<>ToString[\[CapitalUpsilon]\[Tau]z]];
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]z; (* Frequency of mode *)
  \[Omega]1 = 3*\[CapitalUpsilon]\[Tau]*\[Omega]/(2*Sqrt[K0]*\[CapitalUpsilon]t);
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,
                                     SetPrecision[a, OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega], OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                         "precODE" -> OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  R = RCorrection[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],
                         SetPrecision[\[Omega],OptionValue[WorkingPrecision]+5],1,
                         SetPrecision[\[Lambda],OptionValue[WorkingPrecision]+5],
                         SetPrecision[d\[Lambda]d\[Omega],OptionValue[WorkingPrecision]+5],e,p,"precODE"->OptionValue[WorkingPrecision]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  d\[ScriptCapitalC]2d\[Omega] = 4 \[Lambda]^3 d\[Lambda]d\[Omega]+4 \[Lambda]^2 (3 d\[Lambda]d\[Omega] + 10 a (m-2 a \[Omega]))+8 \[Lambda] (d\[Lambda]d\[Omega] (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega])) + 
        48 \[Omega] (a m d\[Lambda]d\[Omega]+6-18 a^3 m \[Omega]+12 a^4 \[Omega]^2+a^2 (d\[Lambda]d\[Omega] \[Omega]+6 m^2));  (* \[Omega] derivative of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  d\[Alpha]d\[Omega] = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*d\[ScriptCapitalC]2d\[Omega] + 256*(2*rplus)^5*((P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2)/\[ScriptCapitalC]2; (* \[Omega] derivative of the constant for horizon fluxes *)
  sumPlus0  = sumPlus1  = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* numbers of steps for wr and wz integration *)
  stepsr = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  stepsz = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[stepsz]<>" steps in wz"];
  zlist = {};(* List for functions of \[Theta] *)
  wr = N[Range[1/2*2Pi/stepsr,(stepsr-1)/2*2Pi/stepsr,2Pi/stepsr],Precision[{a,p,e,x}]];
  rp = r[wr];
  Urp = Outer[Times, Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+K0)], {1,-1,1,-1}];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
  expr = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr]) + n*wr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
  expr1 = Outer[Times, 3*I*\[Omega]/(2*Sqrt[K0])*(\[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr] - orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]r"][wr]), {1,-1,1,-1}];(* Linear part of the exponential term *)
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKdr = 2*rp*\[Omega];
  dKd\[Omega]  = (rp^2+a^2);
  d2Kdrd\[Omega] = 2*rp;
  d\[CapitalDelta] = 2*(rp-1);
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,x}]];
    rp = r[wr];
    Urp = {1,-1,1,-1}*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+K0)];(* Geodesic radial velocity at each quadrant (positive and negative radial and polar velocity) *)
    expr = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with geodesic \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
    expr1 = {1,-1,1,-1}*3*I*\[Omega]/(2*Sqrt[K0])*(\[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr] - orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]r"][wr]);(* Linear part of the exponential term *)
    \[CapitalDelta]  = rp^2-2rp+a^2;
    K  = (rp^2+a^2)*\[Omega]-a*m;
    dKdr = 2*rp*\[Omega];
    dKd\[Omega]  = (rp^2+a^2);
    d2Kdrd\[Omega] = 2*rp;
    d\[CapitalDelta] = 2*(rp-1);
    V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
    dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
    (* In solutions of the radial equation and their r and \[Omega] derivatives *)
    RIn    = R[rp]["R0"]["In"]["R"]; 
    dRIndr   = R[rp]["R0"]["In"]["dR"];
    d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dRInd\[Omega]   = R[rp]["R1"]["In"]["R"]; (* \[Omega]-derivatives of the radial function *)
    d2RIndrd\[Omega]  = R[rp]["R1"]["In"]["dR"];
    d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
    (* Operators D and their r and \[Omega] derivatives *)
    DRIn = dRIndr - I*K/\[CapitalDelta]*RIn;
    dDRInd\[Omega] = d2RIndrd\[Omega] - I*K/\[CapitalDelta]*dRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RIn;
    dDRIndr = d2RIndr2 - I*K/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
    d2DRIndrd\[Omega] = d3RIndr2d\[Omega] - I*K/\[CapitalDelta]*d2RIndrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRInd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
    DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
    dDDRInd\[Omega] = d2DRIndrd\[Omega] - I*K/\[CapitalDelta]*dDRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRIn;
    (* Up solutions of the radial equation and their r and \[Omega] derivatives *)
    RUp    = R[rp]["R0"]["Up"]["R"]; 
    dRUpdr   = R[rp]["R0"]["Up"]["dR"];
    d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
    dRUpd\[Omega]   = R[rp]["R1"]["Up"]["R"]; (* \[Omega]-derivative of the radial function *)
    d2RUpdrd\[Omega]  = R[rp]["R1"]["Up"]["dR"];
    d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
    (* Operators D and their r and \[Omega] derivatives *)
    DRUp = dRUpdr - I*K/\[CapitalDelta]*RUp;
    dDRUpd\[Omega] = d2RUpdrd\[Omega] - I*K/\[CapitalDelta]*dRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RUp;
    dDRUpdr = d2RUpdr2 - I*K/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
    d2DRUpdrd\[Omega] = d3RUpdr2d\[Omega] - I*K/\[CapitalDelta]*d2RUpdrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRUpd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
    DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
    dDDRUpd\[Omega] = d2DRUpdrd\[Omega] - I*K/\[CapitalDelta]*dDRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRUp;
    (* integration over w\[Theta] *)
    For[iz = 1, iz <= stepsz/2, iz++,
      wz = N[(iz-1/2)*2Pi/stepsz,Precision[{a,p,e,x}]];
      If[ir==1,(* functions of only \[Theta] saved to a list in the first step *)
        zp = z[wz];
        \[Theta]p = ArcCos[zp];
        Uzp = {1,1,-1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(K0-a^2*zp^2)];(* Polar geodesic velocity *)
        expz = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][wz])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][wz])+2Pi*k*(iz-1/2)/stepsz)];(* Exponential term with geodesic \[CapitalDelta]tz and \[CapitalDelta]\[Phi]z *)
        expz1 = {1,1,-1,-1}*3*I*\[Omega]/(2*Sqrt[K0])*(\[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][wz] - orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]\[Theta]"][wz]);(* Linear part of the exponential term *)
        sin\[Theta]p = Sqrt[1-zp^2];
        K\[Theta] = a*\[Omega]*(1-zp^2) - m;
        dK\[Theta]d\[Omega] = a*(1-zp^2);
        (* Spin-weighted spheroidal harmonics and their \[Theta]-derivatives *)
        {S,dSd\[Theta]}=SWSH[\[Theta]p]; 
        d2Sd\[Theta]2 = -(zp/sin\[Theta]p)*dSd\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
        (* \[Omega]-derivatives of S *)
        {dSd\[Omega],d2Sd\[Theta]d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
        d3Sd\[Theta]2d\[Omega] = -(zp/sin\[Theta]p)*d2Sd\[Theta]d\[Omega] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega]*(1-zp^2) + 4*a*zp + 2*m*a + d\[Lambda]d\[Omega])*S;
        (* Operators L and their \[Omega]-derivatives *)
        L2S = dSd\[Theta] + (K\[Theta] + 2 zp)*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        dL2Sd\[Theta] = d2Sd\[Theta]2 + (K\[Theta] + 2 zp)*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*S;
        L1L2S = dL2Sd\[Theta] + (K\[Theta] + zp)*L2S/sin\[Theta]p; 
        dL2Sd\[Omega] = d2Sd\[Theta]d\[Omega] + (K\[Theta] + 2 zp)*dSd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
        d2L2Sd\[Theta]d\[Omega] = d3Sd\[Theta]2d\[Omega] + (K\[Theta] + 2 zp)*d2Sd\[Theta]d\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*dSd\[Omega] + (2*a*zp - dK\[Theta]d\[Omega]/sin\[Theta]p^2*zp)*S;
        dL1L2Sd\[Omega] = d2L2Sd\[Theta]d\[Omega] + (K\[Theta] + zp)*dL2Sd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*L2S/sin\[Theta]p; 
        AppendTo[zlist,{zp,Uzp,expz,expz1,sin\[Theta]p,S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]}],
        {zp,Uzp,expz,expz1,sin\[Theta]p,S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]}=zlist[[iz]];
      ];
      (* Integration f functions of r and \[Theta] *)
      \[Zeta] = rp-I*a*zp;
      \[Zeta]bar = rp+I*a*zp;
      vn = -((rp^2+a^2)*En0 - a*Lz0 + Urp)/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
      vmb = (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uzp)/sin\[Theta]p^2;
      (* In solution *)
      {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RIn,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
      {dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega]}   = dFabd\[Omega][\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RIn,DRIn,DDRIn,dRInd\[Omega],dDRInd\[Omega],dDDRInd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
      (* Summing the geodesic ampltude over all quadrants *)
      sumPlus0  += Total[(vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];
      (* Summing the linear part of the ampltude over all quadrants *)
      sumPlus1  += Total[(
          (vn*( GnIn ) + vmb*( GmbIn ))/Sqrt[K0] +
          (vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*(expr1 + expz1) +
          (vn*vn*dFnnInd\[Omega] + vn*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
         )*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}]; 
      (* Up solution *)
      {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RUp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
      {dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega]}   = dFabd\[Omega][\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RUp,DRUp,DDRUp,dRUpd\[Omega],dDRUpd\[Omega],dDDRUpd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
      (* Summing the geodesic ampltude over all quadrants *)
      sumMinus0  += Total[(vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}];
      (* Summing the linear part of the ampltude over all quadrants *)
      sumMinus1  += Total[(
          (vn*GnUp + vmb*GmbUp)/Sqrt[K0] +
          (vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*(expr1 + expz1) +
          (vn*vn*dFnnUpd\[Omega] + vn*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
         )*expr^{1,-1,1,-1}*expz^{1,1,-1,-1}]; 
    ];
  ];
  W = (RIn*dRUpdr - dRIndr*RUp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  W1 = (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/(rp^2-2*rp+a^2)*\[Omega]1; (* Invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalUpsilon]t*W*stepsr*stepsz); (* Amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalUpsilon]t*W*stepsr*stepsz);
  CPlus1  = 2*Pi*(sumPlus1)/(\[CapitalUpsilon]t*W*stepsr*stepsz) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CPlus0; (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1)/(\[CapitalUpsilon]t*W*stepsr*stepsz) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CMinus0;
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->n,
    "\[Omega]"->\[Omega],
    "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"->(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*Abs[CPlus0]^2/(2Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*\[Alpha]*Abs[CMinus0]^2/(2Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"-> (En0*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z) - a*( m - a*\[Omega] ) - 3*\[Omega]*( -\[CapitalUpsilon]\[Tau]z + \[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*\[CapitalUpsilon]tz ))/(2*Sqrt[K0])*Abs[CPlus0]^2/(2Pi*\[Omega]^3)+(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(2Pi*\[Omega]^3),
        "\[ScriptCapitalH]"-> (En0*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z) - a*( m - a*\[Omega] ) - 3*\[Omega]*( -\[CapitalUpsilon]\[Tau]z + \[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*\[CapitalUpsilon]tz ))/(2*Sqrt[K0])*\[Alpha]*Abs[CMinus0]^2/(2Pi*\[Omega]^3)+(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(2Pi*\[Omega]^3)
      |>
    |>,
    "stepsr"->stepsr,
    "steps\[Theta]"->stepsz
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


(* ::Subsection::Closed:: *)
(*Spherical*)


TeukolskySpinModeSpherical[l_?IntegerQ,m_?IntegerQ,k_?IntegerQ,orbitCorrection_,s_]:=Module[{h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],correction,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,
    sumPlus,sumMinus,steps\[Theta],correctionp,i\[Theta],w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,RUprp,dRUprp,ddRUprp,dddRUprp,
    \[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus,CMinus},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_] := 1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+
               a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  If[!PossibleZeroQ[e],Return[$Failed]];
  x = orbitCorrection["Inclination"];
  En = orbitCorrection["En0"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lz0"];
  Kc = orbitCorrection["K0"];
  En1 = orbitCorrection["\[Delta]En"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["\[Delta]Lz"];
  {\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]+s*orbitCorrection["BLFrequenciesCorrection"]; (* BL frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];(* Geodesic coordinate z=cos(\[Theta]) *)
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]]; 
  \[Omega] = m*\[CapitalOmega]\[Phi] + k*\[CapitalOmega]\[Theta]; (* Frequency *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]]; (* Polar and radial functions and the eigenvalue *)
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega]];
  \[Lambda] = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalEigenvalue[-2,l,m,a*\[Omega]];
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2);  (*  TS constant *)
  rplus = 1+Sqrt[1-a^2];  (* Outer horizon *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  sumPlus  = 0; (* Results of the integration are stored in these variables *)
  sumMinus = 0;
  (* number of steps for w\[Theta] integration *)
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  (*Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];*)
  correctionp = Table[correction[N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]]],{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points *)
  rp = p;
  Ur = 0;(* Geodesic radial velocity *)
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  d\[CapitalDelta] = 2*(rp-1);
  dK = 2rp*\[Omega];
  V  = -(K^2+4I*(rp-1)*K)/\[CapitalDelta]+8I*\[Omega]*rp+\[Lambda]; (* Potential in radial Teukolsky equation *)
  dV = -((2K*dK+4I*K+4I*(rp-1)*dK)*\[CapitalDelta]-(K^2+4I*(rp-1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2+8I*\[Omega]; (* derivative of potential in radial Teukolsky equation *)
  RInrp    = R["In"][rp]; (* radial function *)
  dRInrp   = R["In"]'[rp];
  ddRInrp  = (V*RInrp+d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  dddRInrp = 1/\[CapitalDelta] (dV*RInrp+(V+2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
  RUprp    = R["Up"][rp];
  dRUprp   = R["Up"]'[rp];
  ddRUprp  = (V*RUprp+d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
  dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
  For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++,(* integration over w\[Theta] *)
    w\[Theta] = N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]];
    zp = z[w\[Theta]];
    Uz = {1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
    exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
    sin\[Theta]p = Sqrt[1-zp^2];
    S = SWSH[ArcCos[zp],0];  (*  Polar function S(\[Theta](z))  *)
    dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]);
    d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta];(* second derivative of S from Teukolsky equation *)
    d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
             -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
    L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
    dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
    L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
    dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
               1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
    \[Zeta] = rp-I*a*zp;
    \[Zeta]bar = rp+I*a*zp;
    \[CapitalSigma] = rp^2+a^2*zp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
    {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
    vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
    vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
    vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
    vmb = Conjugate[vm];
    Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
    Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
    Snmb = Conjugate[Snm];
    Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
    Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
    Amnn   = vn^2;
    Amnmb  = vn*vmb;
    Ammbmb = vmb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Adnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Admbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Adnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
      -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
    St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
    S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
    rp1 = {correctionp[[i\[Theta]]]["\[Delta]r"],correctionp[[-i\[Theta]]]["\[Delta]r"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
    zp1 = {correctionp[[i\[Theta]]]["\[Delta]z"],-correctionp[[-i\[Theta]]]["\[Delta]z"]};
    Urp1 = {correctionp[[i\[Theta]]]["\[Delta]Ur"],correctionp[[-i\[Theta]]]["\[Delta]Ur"]};
    Uzp1 = {correctionp[[i\[Theta]]]["\[Delta]Uz"],-correctionp[[-i\[Theta]]]["\[Delta]Uz"]};
    \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
    exp1 = {I*(\[Omega]*correctionp[[i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"]),
            I*(\[Omega]*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"])};(* Linear parts of the exponential terms with \[CapitalDelta]t and \[CapitalDelta]\[Phi] *)
    vn1  = -( ((2*rp*rp1)*En) - ((rp^2+a^2)*En - a*Lz + Ur)/(\[CapitalSigma])*\[CapitalSigma]1 + 
             ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]) + Urp1))/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
    vmb1 = ( (-I*(-2*a*zp*zp1*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp])) + Uzp1))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Adnn;
    Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Adnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Admbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
    AmbmbrS  = (Srmb + rp1*vmb)*vmb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
    sumPlus  += Total[s*\[CapitalSigma]*((
      ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
      ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
      ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
      ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
      ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
      ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
      \[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RInrp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRInrp + (Ammbmb*fmbmb2)*ddRInrp)*exp\[Theta]^{1,-1})]; 
    sumMinus += Total[s*\[CapitalSigma]*((
      ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
      ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
      ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
      ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
      ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
      ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
      \[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RUprp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRUprp + (Ammbmb*fmbmb2)*ddRUprp)*exp\[Theta]^{1,-1})];
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  CPlus  = 2*Pi*(1-s*\[CapitalGamma]1/\[CapitalGamma])*sumPlus/(\[CapitalGamma]*W*steps\[Theta]); (* Amplitudes *)
  CMinus = 2*Pi*(1-s*\[CapitalGamma]1/\[CapitalGamma])*sumMinus/(\[CapitalGamma]*W*steps\[Theta]);
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->0,
    "\[Omega]"->\[Omega],
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus,
      "\[ScriptCapitalH]"->CMinus
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "steps\[Theta]"->steps\[Theta]
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


TeukolskySpinModeSphericalCorrectionNum[l_?IntegerQ,m_?IntegerQ,k_?IntegerQ,orbitCorrection_,\[Delta]\[Omega]_]:=Module[{h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,correction,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,
    SWSH,SWSHplus,SWSHminus,R,Rplus,Rminus,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,steps\[Theta],correctionp,i\[Theta],w\[Theta],rp,zp,sin\[Theta]p,
    Ur,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,
    \[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,
    fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_] := 1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+
               a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  x = orbitCorrection["Inclination"];
  En = orbitCorrection["En0"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lz0"];
  Kc = orbitCorrection["K0"];
  En1 = orbitCorrection["\[Delta]En"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["\[Delta]Lz"];
  {\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Geodesic BL frequencies *)
  {\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Coordinate frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];(* Geodesic coordinate z=cos(\[Theta]) *)
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega] = m*\[CapitalOmega]\[Phi] + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + k*\[CapitalOmega]\[Theta]1;
  SWSH = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*\[Omega]];(* Polar and radial functions and the eigenvalue for geodesic frequency and for numerical derivative *)
  SWSHplus = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*(\[Omega]+\[Delta]\[Omega])];
  SWSHminus = SpinWeightedSpheroidalHarmonics`SpinWeightedSpheroidalHarmonicS[-2,l,m,a*(\[Omega]-\[Delta]\[Omega])];
  R = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega]];
  Rplus = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega]+\[Delta]\[Omega]];
  Rminus = Teukolsky`TeukolskyRadial`TeukolskyRadial[-2,l,m,a,\[Omega]-\[Delta]\[Omega]];
  \[Lambda] = R["In"]["Eigenvalue"];
  \[Lambda]1 = (Rplus["In"]["Eigenvalue"] - Rminus["In"]["Eigenvalue"])/(2*\[Delta]\[Omega])*\[Omega]1; (* Linear correction to the eigenvalue *)
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  \[ScriptCapitalC]21 = 4 \[Lambda]^3 \[Lambda]1+4 \[Lambda]^2 (3 \[Lambda]1+10 a (m-2 a \[Omega]) \[Omega]1)+8 \[Lambda] (\[Lambda]1 (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega]) \[Omega]1) + 
        48 \[Omega] (a m \[Lambda]1+6 \[Omega]1-18 a^3 m \[Omega] \[Omega]1+12 a^4 \[Omega]^2 \[Omega]1+a^2 (\[Lambda]1 \[Omega]+6 m^2 \[Omega]1));  (* linear part of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  \[Alpha]1 = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*\[ScriptCapitalC]21 + 256*(2*rplus)^5*(\[Omega]1*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P*\[Omega]1)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P*\[Omega]1)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2*\[Omega]1)/\[ScriptCapitalC]2; (* linear part of the constant for horizon fluxes *)
  sumPlus0 = sumPlus1 = 0; (* Results of the integration are stored in these variables *)
  sumMinus0 = sumMinus1 = 0;
  (* number of steps for w\[Theta] integration *)
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  (*Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];*)
  correctionp = Table[correction[N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]]],{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points *)
  rp = p;
  Ur = 0;(* Geodesic radial velocity *)
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  K1  = (rp^2+a^2)*\[Omega]1;
  d\[CapitalDelta] = 2*(rp-1);
  dK = 2*rp*\[Omega];
  dK1 = 2*rp*\[Omega]1;
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  V1  = -(2*K*K1 + 4I*(rp-1)*K1)/\[CapitalDelta] + 8*I*\[Omega]1*rp + \[Lambda]1; (* Linear part of the potential in radial Teukolsky equation *)
  dV = -((2K*dK + 4I*K + 4I*(rp - 1)*dK)*\[CapitalDelta] - (K^2 + 4I*(rp - 1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2 + 8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
  RInrp    = R["In"][rp]; (* radial function *)
  dRInrp   = R["In"]'[rp];
  ddRInrp  = (V*RInrp + d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  RInrp1   = (Rplus["In"][rp] - Rminus["In"][rp])/(2*\[Delta]\[Omega])*\[Omega]1; (* Linear part of the radial function *)
  dRInrp1  = (Rplus["In"]'[rp] - Rminus["In"]'[rp])/(2*\[Delta]\[Omega])*\[Omega]1;
  ddRInrp1 = (V*RInrp1 + V1*RInrp + d\[CapitalDelta]*dRInrp1)/\[CapitalDelta];
  dddRInrp = 1/\[CapitalDelta] (dV*RInrp + (V + 2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
  RUprp    = R["Up"][rp];
  dRUprp   = R["Up"]'[rp];
  ddRUprp  = (V*RUprp + d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
  RUprp1   = (Rplus["Up"][rp] - Rminus["Up"][rp])/(2*\[Delta]\[Omega])*\[Omega]1;
  dRUprp1  = (Rplus["Up"]'[rp] - Rminus["Up"]'[rp])/(2*\[Delta]\[Omega])*\[Omega]1;
  ddRUprp1 = (V*RUprp1 + V1*RUprp + d\[CapitalDelta]*dRUprp1)/\[CapitalDelta];
  dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
  For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++,(* integration over w\[Theta] *)
    w\[Theta] = N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]];
    zp = z[w\[Theta]];
    Uz = {1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
    exp\[Theta] = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
    sin\[Theta]p = Sqrt[1-zp^2];
    S = SWSH[ArcCos[zp],0];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
    S1 = (SWSHplus[ArcCos[zp],0]-SWSHminus[ArcCos[zp],0])/(2*\[Delta]\[Omega])*\[Omega]1;  (*  Linear part of S(\[Theta](z))  *)
    dSd\[Theta] = (D[SWSH[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]); (* First derivative of S wrt \[Theta] *)
    dS1d\[Theta] = ((D[SWSHplus[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp])-(D[SWSHminus[\[Theta]2,0],\[Theta]2]/.\[Theta]2->ArcCos[zp]))/(2*\[Delta]\[Omega])*\[Omega]1; (* Linear part of first derivative of S wrt \[Theta] *)
    d2Sd\[Theta]2 = -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S-zp/sin\[Theta]p*dSd\[Theta]; (* second derivative of S from Teukolsky equation *)
    d2S1d\[Theta]2 = (-(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda])*S1-(-a^2*2*\[Omega]*\[Omega]1*(1-zp^2)+4*a*zp*\[Omega]1+2*m*a*\[Omega]1+\[Lambda]1)*S-zp/sin\[Theta]p*dS1d\[Theta]);(* Linear part of second derivative of S from Teukolsky equation *)
    d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
             -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
    L2S = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
    L2S1 = dS1d\[Theta]-(S1 (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2) \[Omega]1))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
    dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
    L1L2S = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
    L1L2S1 = d2S1d\[Theta]2+(dS1d\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2) \[Omega]1))/sin\[Theta]p + 
             S1 (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) \[Omega]1 - a^2 (-1+zp^2) 2*\[Omega]*\[Omega]1);  
    dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
               1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
    \[Zeta] = rp-I*a*zp;
    \[Zeta]bar = rp+I*a*zp;
    \[CapitalSigma] = rp^2+a^2*zp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,a,\[Omega]];
    {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
    {fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21} = dfabidS[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dK,K1,dK1,S,L2S,L1L2S,S1,L2S1,L1L2S1,a];
    vn = -((rp^2+a^2)*En - a*Lz + Ur)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
    vl = -((rp^2+a^2)*En - a*Lz - Ur)/(\[CapitalDelta]);
    vm = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
    vmb = Conjugate[vm];
    Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
    Snm  = (\[Zeta]/Sqrt[Kc])*vm*vn;
    Snmb = Conjugate[Snm];
    Slmb = (-(\[Zeta]/Sqrt[Kc]))*vl*vmb;
    Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
    Amnn   = vn^2;
    Amnmb  = vn*vmb;
    Ammbmb = vmb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Adnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Admbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Adnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
      -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
    St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
    S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
    rp1 = {correctionp[[i\[Theta]]]["\[Delta]r"],correctionp[[-i\[Theta]]]["\[Delta]r"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
    zp1 = {correctionp[[i\[Theta]]]["\[Delta]z"],-correctionp[[-i\[Theta]]]["\[Delta]z"]};
    Urp1 = {correctionp[[i\[Theta]]]["\[Delta]Ur"],correctionp[[-i\[Theta]]]["\[Delta]Ur"]};
    Uzp1 = {correctionp[[i\[Theta]]]["\[Delta]Uz"],-correctionp[[-i\[Theta]]]["\[Delta]Uz"]};
    \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
    exp1 = I*{(\[Omega]*correctionp[[ i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[ i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
              (\[Omega]*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};
    vn1  = -( ((2*rp*rp1)*En) - ((rp^2+a^2)*En - a*Lz + Ur)/(\[CapitalSigma])*\[CapitalSigma]1 + 
             ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]) + Urp1))/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
    vmb1 = ( (-I*(-2*a*zp*zp1*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp])) + Uzp1))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Adnn;
    Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Adnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Admbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
    AmbmbrS  = (Srmb + rp1*vmb)*vmb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RInrp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRInrp + Ammbmb*fmbmb2*ddRInrp)*exp\[Theta]^{1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RUprp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRUprp + Ammbmb*fmbmb2*ddRUprp)*exp\[Theta]^{1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*RInrp1 - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*dRInrp1 + Ammbmb*fmbmb2*ddRInrp1
          + (Amnn*fnn01 + Amnmb*fnmb01 + Ammbmb*fmbmb01)*RInrp - (Amnmb*fnmb11 + Ammbmb*fmbmb11)*dRInrp + Ammbmb*fmbmb21*ddRInrp)*exp\[Theta]^{1,-1})]; 
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*RUprp1 - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*dRUprp1 + Ammbmb*fmbmb2*ddRUprp1
          + (Amnn*fnn01 + Amnmb*fnmb01 + Ammbmb*fmbmb01)*RUprp - (Amnmb*fnmb11 + Ammbmb*fmbmb11)*dRUprp + Ammbmb*fmbmb21*ddRUprp)*exp\[Theta]^{1,-1})]; 
  ];
  W = (RInrp*dRUprp - dRInrp*RUprp)/\[CapitalDelta]; (* Invariant Wronskian *)
  W1 = (RInrp1*dRUprp + RInrp*dRUprp1 - dRInrp1*RUprp - dRInrp*RUprp1)/\[CapitalDelta]; (* Linear part of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - W1/W*sumPlus0 )/(\[CapitalGamma]*W*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - W1/W*sumMinus0)/(\[CapitalGamma]*W*steps\[Theta]);
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->0,
    "\[Omega]"->\[Omega],
    "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha],
    "\[Alpha]Correction"->\[Alpha]1,
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(\[Alpha]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "\[Lambda]"->\[Lambda],
    "\[Lambda]Correction"->\[Lambda]1,
    "steps\[Theta]"->steps\[Theta]
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


Options[TeukolskySpinModeSphericalCorrection] = {WorkingPrecision->30};

TeukolskySpinModeSphericalCorrection[l_?IntegerQ,m_?IntegerQ,k_?IntegerQ,orbitCorrection_,orbitDerivatives_,{angparNew_,TeukolskySolverHS1spin_},OptionsPattern[]]:=Module[{
    h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,dEndr,dLzdr,dEndx,dLzdx,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,d\[CapitalOmega]\[Theta]dr,d\[CapitalOmega]\[Phi]dr,d\[CapitalOmega]\[Theta]dx,d\[CapitalOmega]\[Phi]dx,correction,derivatives,z,\[CapitalGamma],\[CapitalGamma]1,d\[CapitalGamma]dr,d\[CapitalGamma]dx,
    \[Omega],\[Omega]prec,\[Omega]1,d\[Omega]dr,d\[Omega]dx,\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega],R,\[ScriptCapitalC]2,d\[ScriptCapitalC]2d\[Omega],rplus,P,\[Epsilon],\[Alpha],d\[Alpha]d\[Omega],
    steps\[Theta],correctionp,derivativesp,sumPlus0,sumMinus0,sumPlus1,sumMinus1,dsumPlusdr,dsumMinusdr,dsumPlusdx,dsumMinusdx,i\[Theta],w\[Theta],
    rp,zp,\[Theta]p,sin\[Theta]p,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dKd\[Omega],dKdr,d2Kdrd\[Omega],V,dVd\[Omega],dVdr,
    RIn,dRIndr,d2RIndr2,d3RIndr3,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,d3RUpdr3,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    \[Theta]2,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],L2S,dL2Sd\[Theta],dL2Sd\[Omega],L1L2S,dL1L2Sd\[Theta],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0d\[Omega],dfnmb0d\[Omega],dfnmb1d\[Omega],dfmbmb0d\[Omega],dfmbmb1d\[Omega],dfmbmb2d\[Omega],
    dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],
    ul,un,um,umb,Sln,Slmb,Snm,Snmb,Smmb,
    Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,
    St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,
    rp1,zp1,dzpdr,d\[Theta]pdr,dzpdx,d\[Theta]pdx,Urp1,Uzp1,dUzpdr,dUzpdx,\[CapitalSigma]1,d\[CapitalSigma]dr,d\[CapitalSigma]dx,exp1,dexpdr,dexpdx,
    un1,dundr,dundx,umb1,dumbdr,dumbdx,
    Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,W,dWd\[Omega],
    CPlus0,CMinus0,CPlus1,CMinus1,dCPlusdr,dCMinusdr,dCPlusdx,dCMinusdx},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_] := 1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+
               a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  (*Print["Calculating l = "<>ToString[l]<>", m = "<>ToString[m]<>", k = "<>ToString[k]<>" mode"];*)
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  If[e!=0,Return[$Failed]];
  x = orbitCorrection["Inclination"];
  En = orbitCorrection["En0"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lz0"];
  Kc = orbitCorrection["K0"];
  En1 = orbitCorrection["\[Delta]En"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["\[Delta]Lz"];
  dEndr = orbitDerivatives["dEndr"]; (* Derivatives of the constants of motion *)
  dLzdr = orbitDerivatives["dLzdr"];
  dEndx = orbitDerivatives["dEndx"]; 
  dLzdx = orbitDerivatives["dLzdx"];
  {\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Coordinate frequencies *)
  {\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the coordinate frequencies *)
  {d\[CapitalOmega]\[Theta]dr,d\[CapitalOmega]\[Phi]dr} = orbitDerivatives["dBLFrequenciesdr"]; (* Derivatives of the coordinate frequencies *)
  {d\[CapitalOmega]\[Theta]dx,d\[CapitalOmega]\[Phi]dx} = orbitDerivatives["dBLFrequenciesdx"];
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  derivatives = orbitDerivatives["OrbitDerivatives"]; (* function containing derivatives of the trajectory *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction and derivatives *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]]; 
  d\[CapitalGamma]dr = orbitDerivatives["dMinoFrequenciesdr"][[1]]; 
  d\[CapitalGamma]dx = orbitDerivatives["dMinoFrequenciesdx"][[1]]; 
  \[Omega] = m*\[CapitalOmega]\[Phi] + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction and derivatives *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  \[Omega]prec = m*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"]
        + k*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)"]; (* Frequency with higher precision *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + k*\[CapitalOmega]\[Theta]1;
  d\[Omega]dr = m*d\[CapitalOmega]\[Phi]dr + k*d\[CapitalOmega]\[Theta]dr;
  d\[Omega]dx = m*d\[CapitalOmega]\[Phi]dx + k*d\[CapitalOmega]\[Theta]dx;
  (*Print["Calculating angular functions"<>ToString@AbsoluteTiming[{\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1];][[1]]];*)(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]} = angparNew[-2,l,m,
                                     SetPrecision[a,    OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega]prec,OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                           "precODE" -> OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)];(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  (*Print["Calculating radial functions"<>ToString@AbsoluteTiming[*)R = TeukolskySolverHS1spin[p,-2,l,m,
                                                                SetPrecision[a,    OptionValue[WorkingPrecision]+5],
                                                                SetPrecision[\[Omega]prec,OptionValue[WorkingPrecision]+5],1,
                                                                SetPrecision[\[Lambda],    OptionValue[WorkingPrecision]+5],
                                                                SetPrecision[d\[Lambda]d\[Omega], OptionValue[WorkingPrecision]+5],
                                                                       "precODE" -> OptionValue[WorkingPrecision]](*][[1]]]*);
  (*R = TeukolskySolverHS1spin[p,-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1,\[Lambda],d\[Lambda]d\[Omega]];*)
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  d\[ScriptCapitalC]2d\[Omega] = 4 \[Lambda]^3 d\[Lambda]d\[Omega]+4 \[Lambda]^2 (3 d\[Lambda]d\[Omega]+10 a (m-2 a \[Omega]))+8 \[Lambda] (d\[Lambda]d\[Omega] (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega])) + 
          48 \[Omega] (a m d\[Lambda]d\[Omega]+6-18 a^3 m \[Omega]+12 a^4 \[Omega]^2+a^2 (d\[Lambda]d\[Omega] \[Omega]+6 m^2));  (* \[Omega]-derivative of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  d\[Alpha]d\[Omega] = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*d\[ScriptCapitalC]2d\[Omega] + 256*(2*rplus)^5*((P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2)/\[ScriptCapitalC]2; (* \[Omega]-derivative of the constant for horizon fluxes *)
  (* number of steps for w\[Theta] integration *)
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/2]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/2]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  (*Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];*)
  correctionp = Table[correction[N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]]],{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points *)
  derivativesp = Table[derivatives[N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]]],{i\[Theta],1,steps\[Theta]/2}];(* derivatives of the trajectory at all points *)
  rp = p;
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKd\[Omega]  = (rp^2+a^2);
  d\[CapitalDelta] = 2*(rp-1);
  dKdr = 2*rp*\[Omega];
  d2Kdrd\[Omega] = 2*rp;
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  dVdr = -((2K*dKdr + 4I*K + 4I*(rp - 1)*dKdr)*\[CapitalDelta] - (K^2 + 4I*(rp - 1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2 + 8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
  {{{RIn,dRIndr},{dRInd\[Omega],d2RIndrd\[Omega]}},{{RUp,dRUpdr},{dRUpd\[Omega],d2RUpdrd\[Omega]}}} = R; (* Radial function and its derivatives *)
  d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from the Teukolsky equation *)
  d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  
  d3RIndr3 = 1/\[CapitalDelta] (dVdr*RIn + (V + 2)*dRIndr);  (* third derivative of radial function from the Teukolsky equation *)
  d3RUpdr3 = 1/\[CapitalDelta] (dVdr*RUp + (V + 2)*dRUpdr);
  d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
  d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
  {sumPlus0, sumMinus0, sumPlus1, sumMinus1, dsumPlusdr, dsumMinusdr, dsumPlusdx, dsumMinusdx} = Sum[(* integration over w\[Theta] *)
    w\[Theta]     = N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]];
    zp     = z[w\[Theta]];
    Uz     = {1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
    exp\[Theta]   = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
    sin\[Theta]p  = Sqrt[1-zp^2];
    \[Theta]p     = ArcCos[zp];
    {S,dSd\[Theta](*,d2Sd\[Theta]2*)}=SWSH[\[Theta]p];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
    d2Sd\[Theta]2 = -(zp/sin\[Theta]p)*dSd\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
    {dSd\[Omega],d2Sd\[Theta]d\[Omega](*,d3Sd\[Theta]2d\[Omega]*)}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
    d3Sd\[Theta]2d\[Omega] = -(zp/sin\[Theta]p)*d2Sd\[Theta]d\[Omega] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega]*(1-zp^2) + 4*a*zp + 2*m*a + d\[Lambda]d\[Omega])*S;
    d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
             -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from the TE *)
    L2S    = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p; (* Operators acting on S(\[Theta]) *)
    L1L2S  = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
    dL2Sd\[Omega] = (d2Sd\[Theta]d\[Omega]-(dSd\[Omega] (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2)))/sin\[Theta]p); (* \[Theta] and \[Omega]-derivatives of the operator *)
    dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
    dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
               1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
    dL1L2Sd\[Omega] = (d3Sd\[Theta]2d\[Omega]+(d2Sd\[Theta]d\[Omega] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2)))/sin\[Theta]p + 
               dSd\[Omega] (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) - a^2 (-1+zp^2) 2*\[Omega]));  
    \[Zeta] = rp-I*a*zp;
    \[Zeta]bar = rp+I*a*zp;
    \[CapitalSigma] = rp^2+a^2*zp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,a,\[Omega]];
    {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
    {dfnn0d\[Omega],dfnmb0d\[Omega],dfnmb1d\[Omega],dfmbmb0d\[Omega],dfmbmb1d\[Omega],dfmbmb2d\[Omega]} = dfabid\[Omega][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega],a];
    un = -((rp^2+a^2)*En - a*Lz)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
    ul = -((rp^2+a^2)*En - a*Lz)/(\[CapitalDelta]);
    um = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
    umb = Conjugate[um];
    (* Dipole term *)
    Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
    Snm  = (\[Zeta]/Sqrt[Kc])*um*un;
    Snmb = Conjugate[Snm];
    Slmb = (-(\[Zeta]/Sqrt[Kc]))*ul*umb;
    Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
    Amnn   = un^2;
    Amnmb  = un*umb;
    Ammbmb = umb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Adnn  = (Scd\[Gamma]ndc*un-Sln*2*Re[gamma]*un-2*Re[Snmb*((Conjugate[alpha]+beta)*un-mu*um)]);
    Admbmb= (Scd\[Gamma]mbdc*umb-Snmb*(-pi*ul)-Slmb*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)+Smmb*(-(-alpha+Conjugate[beta])*umb));
    Adnmb = (Scd\[Gamma]ndc*umb+Scd\[Gamma]mbdc*un-Sln*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)-Snmb*(Conjugate[rho]*un-mu*ul-(Conjugate[alpha]-beta)*umb)
      -Snm*(-(-alpha+Conjugate[beta])*umb)-Snmb*(-Conjugate[pi]*umb-pi*um)-Slmb*(2*Re[gamma]*un)+Smmb*((alpha+Conjugate[beta])*un-Conjugate[mu]*umb))/2;
    St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
    S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
    (* Derivatives of the trajectory *)
    rp1 = {correctionp[[i\[Theta]]]["\[Delta]r"], correctionp[[-i\[Theta]]]["\[Delta]r"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
    zp1 = {correctionp[[i\[Theta]]]["\[Delta]z"],-correctionp[[-i\[Theta]]]["\[Delta]z"]};
    Urp1 = {correctionp[[i\[Theta]]]["\[Delta]Ur"], correctionp[[-i\[Theta]]]["\[Delta]Ur"]};
    Uzp1 = {correctionp[[i\[Theta]]]["\[Delta]Uz"],-correctionp[[-i\[Theta]]]["\[Delta]Uz"]};
    dzpdr = {derivativesp[[i\[Theta]]]["dzdr"],-derivativesp[[-i\[Theta]]]["dzdr"]}; (* Derivatives of the coordinates and four-velocity for each quadrant *)
    dzpdx = {derivativesp[[i\[Theta]]]["dzdx"],-derivativesp[[-i\[Theta]]]["dzdx"]};
    d\[Theta]pdr = -dzpdr/sin\[Theta]p;
    d\[Theta]pdx = -dzpdx/sin\[Theta]p;
    dUzpdr = {derivativesp[[i\[Theta]]]["dUzdr"],-derivativesp[[-i\[Theta]]]["dUzdr"]};
    dUzpdx = {derivativesp[[i\[Theta]]]["dUzdx"],-derivativesp[[-i\[Theta]]]["dUzdx"]};
    \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction and derivatives of \[CapitalSigma] *)
    d\[CapitalSigma]dr = 2*(rp + a^2*zp*dzpdr);
    d\[CapitalSigma]dx = 2*a^2*zp*dzpdx;
    exp1 = I*{(\[Omega]*correctionp[[ i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[ i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
              (\[Omega]*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]}; (* Linear corrections and derivatives of the exponantial term *)
    dexpdr = I*{(\[Omega]*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]tdr"]-m*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]\[Phi]dr"])+d\[Omega]dr*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
                (\[Omega]*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]tdr"]-m*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]\[Phi]dr"])-d\[Omega]dr*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};
    dexpdx = I*{(\[Omega]*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]tdx"]-m*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]\[Phi]dx"])+d\[Omega]dx*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
                (\[Omega]*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]tdx"]-m*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]\[Phi]dx"])-d\[Omega]dx*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};
    un1  = -( ((2*rp*rp1)*En) - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*\[CapitalSigma]1 + 
             ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]) + Urp1))/(2*\[CapitalSigma]);(* Linear parts and derivatives of the four-velocity in Kinnersley tetrad *)
    dundr  = -( ((2*rp)*En) - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*d\[CapitalSigma]dr + 
             ((rp^2+a^2)*dEndr - a*dLzdr))/(2*\[CapitalSigma]);
    dundx  = -( - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*d\[CapitalSigma]dx + 
             ((rp^2+a^2)*dEndx - a*dLzdx))/(2*\[CapitalSigma]);
    umb1 = ( (-I*(-2*a*zp*zp1*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp])) + Uzp1))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    dumbdr = ( (-I*(-2*a*zp*dzpdr*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*dzpdr/sin\[Theta]p^2 + (1-I*a*dzpdr)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*dEndr - dLzdr) + dUzpdr))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    dumbdx = ( (-I*(-2*a*zp*dzpdx*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*dzpdx/sin\[Theta]p^2 + (-I*a*dzpdx)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*dEndx - dLzdx) + dUzpdx))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*un + 2*un1)*un + Adnn; (* Source term *)
    Annt\[Phi]S = (St\[Phi]n + exp1*un)*un;
    AnnrS  = (Srn + rp1*un)*un;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*un)*un;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*un*umb + un1*umb + un*umb1 + Adnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*umb + St\[Phi]mb*un)/2 + exp1*un*umb);
    AnmbrS  = ((Srn*umb + Srmb*un)/2 + rp1*un*umb);
    Anmb\[Theta]S  = ((S\[Theta]n*umb + S\[Theta]mb*un)/2 - zp1/sin\[Theta]p*un*umb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*umb + 2*umb1)*umb + Admbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*umb)*umb;
    AmbmbrS  = (Srmb + rp1*umb)*umb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*umb)*umb;
      {Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RIn - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRIndr + Ammbmb*fmbmb2*d2RIndr2)*exp\[Theta]^{1,-1})], (* Totral of all quadrants *)
       Total[\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RUp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRUpdr + Ammbmb*fmbmb2*d2RUpdr2)*exp\[Theta]^{1,-1})], 
       Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RIn*fnn0 + AnnrS*(dRIndr*fnn0 + RIn*dfnn0dr) + Ann\[Theta]S*RIn*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RIn*fnmb0 + AnmbrS*( RIn*dfnmb0dr +  dRIndr*fnmb0) + Anmb\[Theta]S* RIn*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRIndr*fnmb1 + AnmbrS*(dRIndr*dfnmb1dr + d2RIndr2*fnmb1) + Anmb\[Theta]S*dRIndr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RIn*fmbmb0 + AmbmbrS*(  dRIndr*fmbmb0 +   RIn*dfmbmb0dr) + Ambmb\[Theta]S*  RIn*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRIndr*fmbmb1 + AmbmbrS*( d2RIndr2*fmbmb1 +  dRIndr*dfmbmb1dr) + Ambmb\[Theta]S* dRIndr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RIndr2*fmbmb2 + AmbmbrS*(d3RIndr3*fmbmb2 + d2RIndr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RIndr2*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*dRInd\[Omega] - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*d2RIndrd\[Omega] + Ammbmb*fmbmb2*d3RIndr2d\[Omega]
          + (Amnn*dfnn0d\[Omega] + Amnmb*dfnmb0d\[Omega] + Ammbmb*dfmbmb0d\[Omega])*RIn - (Amnmb*dfnmb1d\[Omega] + Ammbmb*dfmbmb1d\[Omega])*dRIndr + Ammbmb*dfmbmb2d\[Omega]*d2RIndr2)*\[Omega]1*exp\[Theta]^{1,-1})], 
       Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUp*fnn0 + AnnrS*(dRUpdr*fnn0 + RUp*dfnn0dr) + Ann\[Theta]S*RUp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUp*fnmb0 + AnmbrS*( RUp*dfnmb0dr +  dRUpdr*fnmb0) + Anmb\[Theta]S* RUp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUpdr*fnmb1 + AnmbrS*(dRUpdr*dfnmb1dr + d2RUpdr2*fnmb1) + Anmb\[Theta]S*dRUpdr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUp*fmbmb0 + AmbmbrS*(  dRUpdr*fmbmb0 +   RUp*dfmbmb0dr) + Ambmb\[Theta]S*  RUp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUpdr*fmbmb1 + AmbmbrS*( d2RUpdr2*fmbmb1 +  dRUpdr*dfmbmb1dr) + Ambmb\[Theta]S* dRUpdr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RUpdr2*fmbmb2 + AmbmbrS*(d3RUpdr3*fmbmb2 + d2RUpdr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RUpdr2*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*dRUpd\[Omega] - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*d2RUpdrd\[Omega] + Ammbmb*fmbmb2*d3RUpdr2d\[Omega]
          + (Amnn*dfnn0d\[Omega] + Amnmb*dfnmb0d\[Omega] + Ammbmb*dfmbmb0d\[Omega])*RUp - (Amnmb*dfnmb1d\[Omega] + Ammbmb*dfmbmb1d\[Omega])*dRUpdr + Ammbmb*dfmbmb2d\[Omega]*d2RUpdr2)*\[Omega]1*exp\[Theta]^{1,-1})], 
      Total[\[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Amnn*fnn0*RIn + Amnmb*(fnmb0*RIn - fnmb1*dRIndr) + Ammbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundr*fnn0*RIn + (un*dumbdr+dundr*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdr*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Amnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RIn + fnn0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)) + 
               Amnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RIn + fnmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRIndr - fnmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr)) + 
              Ammbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RIn + fmbmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRIndr - fmbmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RIndr2 + fmbmb2*(d3RIndr3 + d3RIndr2d\[Omega]*d\[Omega]dr)))) 
              )*exp\[Theta]^{1,-1}], 
       Total[\[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Amnn*fnn0*RUp + Amnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Ammbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundr*fnn0*RUp + (un*dumbdr+dundr*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdr*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Amnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RUp + fnn0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)) + 
               Amnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RUp + fnmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fnmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr)) + 
              Ammbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RUp + fmbmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fmbmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RUpdr2 + fmbmb2*(d3RUpdr3 + d3RUpdr2d\[Omega]*d\[Omega]dr)))) 
              )*exp\[Theta]^{1,-1}], 
      Total[\[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Amnn*fnn0*RIn + Amnmb*(fnmb0*RIn - fnmb1*dRIndr) + Ammbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundx*fnn0*RIn + (un*dumbdx+dundx*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdx*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Amnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RIn + fnn0*(dRInd\[Omega]*d\[Omega]dx)) + 
               Amnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RIn + fnmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRIndr - fnmb1*(d2RIndrd\[Omega]*d\[Omega]dx)) + 
              Ammbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RIn + fmbmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfmbmb1d\[Theta]*d\[Theta]pdx + dfmbmb1d\[Omega]*d\[Omega]dx)*dRIndr - fmbmb1*(d2RIndrd\[Omega]*d\[Omega]dx) + 
                      (dfmbmb2d\[Theta]*d\[Theta]pdx + dfmbmb2d\[Omega]*d\[Omega]dx)*d2RIndr2 + fmbmb2*(d3RIndr2d\[Omega]*d\[Omega]dx)))) 
              )*exp\[Theta]^{1,-1}], 
       Total[\[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Amnn*fnn0*RUp + Amnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Ammbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundx*fnn0*RUp + (un*dumbdx+dundx*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdx*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Amnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RUp + fnn0*(dRUpd\[Omega]*d\[Omega]dx)) + 
               Amnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RUp + fnmb0*(dRUpd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRUpdr - fnmb1*(d2RUpdrd\[Omega]*d\[Omega]dx)) + 
              Ammbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RUp + fmbmb0*(dRUpd\[Omega]*d\[Omega]dx)
                    - (dfmbmb1d\[Theta]*d\[Theta]pdx + dfmbmb1d\[Omega]*d\[Omega]dx)*dRUpdr - fmbmb1*(d2RUpdrd\[Omega]*d\[Omega]dx) + 
                      (dfmbmb2d\[Theta]*d\[Theta]pdx + dfmbmb2d\[Omega]*d\[Omega]dx)*d2RUpdr2 + fmbmb2*(d3RUpdr2d\[Omega]*d\[Omega]dx)))) 
              )*exp\[Theta]^{1,-1}]
       },
    {i\[Theta], 1, steps\[Theta]/2}
  ];
  W = (RIn*dRUpdr - dRIndr*RUp)/\[CapitalDelta]; (* Invariant Wronskian *)
  dWd\[Omega]= (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/\[CapitalDelta]; (* \[Omega]-derivative of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*\[Omega]1/W*sumPlus0 )/(\[CapitalGamma]*W*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*\[Omega]1/W*sumMinus0)/(\[CapitalGamma]*W*steps\[Theta]);
  dCPlusdr  = 2*Pi*(dsumPlusdr  - d\[CapitalGamma]dr/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*d\[Omega]dr/W*sumPlus0 )/(\[CapitalGamma]*W*steps\[Theta]); (* Derivatives of the amplitudes *)
  dCMinusdr = 2*Pi*(dsumMinusdr - d\[CapitalGamma]dr/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*d\[Omega]dr/W*sumMinus0)/(\[CapitalGamma]*W*steps\[Theta]);
  dCPlusdx  = 2*Pi*(dsumPlusdx  - d\[CapitalGamma]dx/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*d\[Omega]dx/W*sumPlus0 )/(\[CapitalGamma]*W*steps\[Theta]);
  dCMinusdx = 2*Pi*(dsumMinusdx - d\[CapitalGamma]dx/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*d\[Omega]dx/W*sumMinus0)/(\[CapitalGamma]*W*steps\[Theta]);
  (*Print[NumberForm[Re[CPlus0],16]];*)
  <|
    "l" -> l,
    "m" -> m,
    "k" -> k,
    "n" -> 0,
    "\[Omega]" -> \[Omega],
    "\[Omega]Correction" -> \[Omega]1,
    "\[Omega]Derivatives" -> <|
      "r" -> d\[Omega]dr,
      "x" -> d\[Omega]dx
    |>,
    "Amplitudes" -> <|
      "\[ScriptCapitalI]" -> CPlus0,
      "\[ScriptCapitalH]" -> CMinus0
    |>,
    "AmplitudesCorrection" -> <|
      "\[ScriptCapitalI]" -> CPlus1,
      "\[ScriptCapitalH]" -> CMinus1
    |>,
    "AmplitudesDerivatives" -> <|
      "r" -> <|
        "\[ScriptCapitalI]" -> dCPlusdr,
        "\[ScriptCapitalH]" -> dCMinusdr
      |>,
      "x" -> <|
        "\[ScriptCapitalI]" -> dCPlusdx,
        "\[ScriptCapitalH]" -> dCMinusdx
      |>
    |>,
    "\[Alpha]" -> \[Alpha],
    "d\[Alpha]d\[Omega]" -> d\[Alpha]d\[Omega],
    "S" -> N[SWSH[Pi/2,0][[1]]],
    "Fluxes" -> <|
      "Energy" -> <|
        "\[ScriptCapitalI]" -> Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]" -> \[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum" -> <|
        "\[ScriptCapitalI]" -> Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]" -> \[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection" -> <|
      "Energy" -> <|
        "\[ScriptCapitalI]" -> (2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]" -> \[Alpha]*(2*Re[CMinus1*Conjugate[CMinus0]] + (d\[Alpha]d\[Omega]/\[Alpha]*Abs[CMinus0]^2 - 2*Abs[CMinus0]^2/\[Omega])*\[Omega]1)/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum" -> <|
        "\[ScriptCapitalI]" -> (2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesDerivatives" -> <|
      "r" -> <|
        "Energy"  ->  <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdr*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*d\[Omega]dr/\[Omega])/(4Pi*\[Omega]^2),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dr/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdr*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*d\[Omega]dr/\[Omega])/(4Pi*\[Omega]^2)
        |>,
        "AngularMomentum" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdr*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*d\[Omega]dr/\[Omega])*m/(4Pi*\[Omega]^3),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dr/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdr*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*d\[Omega]dr/\[Omega])*m/(4Pi*\[Omega]^3)
        |>
      |>,
      "x" -> <|
        "Energy" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdx*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*d\[Omega]dx/\[Omega])/(4Pi*\[Omega]^2),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dx/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdx*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*d\[Omega]dx/\[Omega])/(4Pi*\[Omega]^2)
        |>,
        "AngularMomentum" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdx*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*d\[Omega]dx/\[Omega])*m/(4Pi*\[Omega]^3),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dx/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdx*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*d\[Omega]dx/\[Omega])*m/(4Pi*\[Omega]^3)
        |>
      |>
    |>,
    "\[Lambda]" -> \[Lambda],
    "d\[Lambda]d\[Omega]" -> d\[Lambda]d\[Omega],
    "steps\[Theta]" -> steps\[Theta]
  |>
]


Options[TeukolskySpinModeSphericalCorrectionNew] = {WorkingPrecision->30};

TeukolskySpinModeSphericalCorrectionNew[l_?IntegerQ,m_?IntegerQ,k_?IntegerQ,orbitCorrection_,orbitDerivatives_,{angparNew_,TeukolskySolverHS1spin_},OptionsPattern[]]:=Module[{
    h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,dEndr,dLzdr,dEndx,dLzdx,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,d\[CapitalOmega]\[Theta]dr,d\[CapitalOmega]\[Phi]dr,d\[CapitalOmega]\[Theta]dx,d\[CapitalOmega]\[Phi]dx,correction,derivatives,z,\[CapitalGamma],\[CapitalGamma]1,d\[CapitalGamma]dr,d\[CapitalGamma]dx,
    \[Omega],\[Omega]prec,\[Omega]1,d\[Omega]dr,d\[Omega]dx,\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega],R,\[ScriptCapitalC]2,d\[ScriptCapitalC]2d\[Omega],rplus,P,\[Epsilon],\[Alpha],d\[Alpha]d\[Omega],
    steps\[Theta],correctionp,derivativesp,sumPlus0,sumMinus0,sumPlus1,sumMinus1,dsumPlusdr,dsumMinusdr,dsumPlusdx,dsumMinusdx,i\[Theta],w\[Theta],
    rp,zp,\[Theta]p,sin\[Theta]p,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dKd\[Omega],dKdr,d2Kdrd\[Omega],V,dVd\[Omega],dVdr,
    RIn,dRIndr,d2RIndr2,d3RIndr3,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,d3RUpdr3,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    \[Theta]2,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],L2S,dL2Sd\[Theta],dL2Sd\[Omega],L1L2S,dL1L2Sd\[Theta],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0d\[Omega],dfnmb0d\[Omega],dfnmb1d\[Omega],dfmbmb0d\[Omega],dfmbmb1d\[Omega],dfmbmb2d\[Omega],
    dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],
    ul,un,um,umb,Sln,Slmb,Snm,Snmb,Smmb,
    Amnn,Amnmb,Ammbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Adnn,Adnmb,Admbmb,
    St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,
    rp1,zp1,dzpdr,d\[Theta]pdr,dzpdx,d\[Theta]pdx,Urp1,Uzp1,dUzpdr,dUzpdx,\[CapitalSigma]1,d\[CapitalSigma]dr,d\[CapitalSigma]dx,exp1,dexpdr,dexpdx,
    un1,dundr,dundx,umb1,dumbdr,dumbdx,
    Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,W,dWd\[Omega],
    CPlus0,CMinus0,CPlus1,CMinus1,dCPlusdr,dCMinusdr,dCPlusdx,dCMinusdx,u\[Theta],dw\[Theta]du\[Theta]},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_] := 1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+
               a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  (*Print["Calculating l = "<>ToString[l]<>", m = "<>ToString[m]<>", k = "<>ToString[k]<>" mode"];*)
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  If[e!=0,Return[$Failed]];
  x = orbitCorrection["Inclination"];
  En = orbitCorrection["En0"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lz0"];
  Kc = orbitCorrection["K0"];
  En1 = orbitCorrection["\[Delta]En"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["\[Delta]Lz"];
  dEndr = orbitDerivatives["dEndr"]; (* Derivatives of the constants of motion *)
  dLzdr = orbitDerivatives["dLzdr"];
  dEndx = orbitDerivatives["dEndx"]; 
  dLzdx = orbitDerivatives["dLzdx"];
  {\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Coordinate frequencies *)
  {\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the coordinate frequencies *)
  {d\[CapitalOmega]\[Theta]dr,d\[CapitalOmega]\[Phi]dr} = orbitDerivatives["dBLFrequenciesdr"]; (* Derivatives of the coordinate frequencies *)
  {d\[CapitalOmega]\[Theta]dx,d\[CapitalOmega]\[Phi]dx} = orbitDerivatives["dBLFrequenciesdx"];
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  derivatives = orbitDerivatives["OrbitDerivatives"]; (* function containing derivatives of the trajectory *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction and derivatives *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]]; 
  d\[CapitalGamma]dr = orbitDerivatives["dMinoFrequenciesdr"][[1]]; 
  d\[CapitalGamma]dx = orbitDerivatives["dMinoFrequenciesdx"][[1]]; 
  \[Omega] = m*\[CapitalOmega]\[Phi] + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction and derivatives *)
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  \[Omega]prec = m*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"]
        + k*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)"]; (* Frequency with higher precision *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + k*\[CapitalOmega]\[Theta]1;
  d\[Omega]dr = m*d\[CapitalOmega]\[Phi]dr + k*d\[CapitalOmega]\[Theta]dr;
  d\[Omega]dx = m*d\[CapitalOmega]\[Phi]dx + k*d\[CapitalOmega]\[Theta]dx;
  (*Print["Calculating angular functions"<>ToString@AbsoluteTiming[{\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1];][[1]]];*)(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]} = angparNew[-2,l,m,
                                     SetPrecision[a,    OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega]prec,OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                           "precODE" -> OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)];(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  (*Print["Calculating radial functions"<>ToString@AbsoluteTiming[*)R = TeukolskySolverHS1spin[p,-2,l,m,
                                                                SetPrecision[a,    OptionValue[WorkingPrecision]+5],
                                                                SetPrecision[\[Omega]prec,OptionValue[WorkingPrecision]+5],1,
                                                                SetPrecision[\[Lambda],    OptionValue[WorkingPrecision]+5],
                                                                SetPrecision[d\[Lambda]d\[Omega], OptionValue[WorkingPrecision]+5],
                                                                       "precODE" -> OptionValue[WorkingPrecision]](*][[1]]]*);
  (*R = TeukolskySolverHS1spin[p,-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1,\[Lambda],d\[Lambda]d\[Omega]];*)
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  d\[ScriptCapitalC]2d\[Omega] = 4 \[Lambda]^3 d\[Lambda]d\[Omega]+4 \[Lambda]^2 (3 d\[Lambda]d\[Omega]+10 a (m-2 a \[Omega]))+8 \[Lambda] (d\[Lambda]d\[Omega] (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega])) + 
          48 \[Omega] (a m d\[Lambda]d\[Omega]+6-18 a^3 m \[Omega]+12 a^4 \[Omega]^2+a^2 (d\[Lambda]d\[Omega] \[Omega]+6 m^2));  (* \[Omega]-derivative of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  d\[Alpha]d\[Omega] = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*d\[ScriptCapitalC]2d\[Omega] + 256*(2*rplus)^5*((P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2)/\[ScriptCapitalC]2; (* \[Omega]-derivative of the constant for horizon fluxes *)
  (* number of steps for w\[Theta] integration *)
  steps\[Theta] = Max[32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/2]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/2]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  (*Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];*)
  correctionp = Table[correction[N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]]],{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points *)
  derivativesp = Table[derivatives[N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]]],{i\[Theta],1,steps\[Theta]/2}];(* derivatives of the trajectory at all points *)
  rp = p;
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKd\[Omega]  = (rp^2+a^2);
  d\[CapitalDelta] = 2*(rp-1);
  dKdr = 2*rp*\[Omega];
  d2Kdrd\[Omega] = 2*rp;
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  dVdr = -((2K*dKdr + 4I*K + 4I*(rp - 1)*dKdr)*\[CapitalDelta] - (K^2 + 4I*(rp - 1)*K)*d\[CapitalDelta])/\[CapitalDelta]^2 + 8I*\[Omega]; (* derivative of potential in radial Teukolsky equation wrt r *)
  {{{RIn,dRIndr},{dRInd\[Omega],d2RIndrd\[Omega]}},{{RUp,dRUpdr},{dRUpd\[Omega],d2RUpdrd\[Omega]}}} = R; (* Radial function and its derivatives *)
  d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from the Teukolsky equation *)
  d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  
  d3RIndr3 = 1/\[CapitalDelta] (dVdr*RIn + (V + 2)*dRIndr);  (* third derivative of radial function from the Teukolsky equation *)
  d3RUpdr3 = 1/\[CapitalDelta] (dVdr*RUp + (V + 2)*dRUpdr);
  d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
  d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
  {sumPlus0, sumMinus0, sumPlus1, sumMinus1, dsumPlusdr, dsumMinusdr, dsumPlusdx, dsumMinusdx} = Sum[(* integration over w\[Theta] *)
    u\[Theta]     = N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]];
    w\[Theta]     = u\[Theta] - Sqrt[1-x^2]*Sin[2*u\[Theta]]/2;
    dw\[Theta]du\[Theta] = 1 - Sqrt[1-x^2]*Cos[2*u\[Theta]];
    zp     = z[w\[Theta]];
    Uz     = {1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
    exp\[Theta]   = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+k*w\[Theta])];
    sin\[Theta]p  = Sqrt[1-zp^2];
    \[Theta]p     = ArcCos[zp];
    {S,dSd\[Theta](*,d2Sd\[Theta]2*)}=SWSH[\[Theta]p];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
    d2Sd\[Theta]2 = -(zp/sin\[Theta]p)*dSd\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
    {dSd\[Omega],d2Sd\[Theta]d\[Omega](*,d3Sd\[Theta]2d\[Omega]*)}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
    d3Sd\[Theta]2d\[Omega] = -(zp/sin\[Theta]p)*d2Sd\[Theta]d\[Omega] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega]*(1-zp^2) + 4*a*zp + 2*m*a + d\[Lambda]d\[Omega])*S;
    d3Sd\[Theta]3 = -(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
             -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from the TE *)
    L2S    = dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p; (* Operators acting on S(\[Theta]) *)
    L1L2S  = d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
    dL2Sd\[Omega] = (d2Sd\[Theta]d\[Omega]-(dSd\[Omega] (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2)))/sin\[Theta]p); (* \[Theta] and \[Omega]-derivatives of the operator *)
    dL2Sd\[Theta] = d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
    dL1L2Sd\[Theta] = d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
               1/sin\[Theta]p^3 (d2Sd\[Theta]2 (-1+zp^2) (2 m-3 zp+2 a (-1+zp^2) \[Omega])+2 S (-m^2 zp+m (1+zp^2)+a (-1+zp^2)^2 \[Omega] (-2+a zp \[Omega])));
    dL1L2Sd\[Omega] = (d3Sd\[Theta]2d\[Omega]+(d2Sd\[Theta]d\[Omega] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2)))/sin\[Theta]p + 
               dSd\[Omega] (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) - a^2 (-1+zp^2) 2*\[Omega]));  
    \[Zeta] = rp-I*a*zp;
    \[Zeta]bar = rp+I*a*zp;
    \[CapitalSigma] = rp^2+a^2*zp^2;
    {fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2} = fabi[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,a];
    {dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr} = dfabidr[\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,a,\[Omega]];
    {dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta]} = dfabid\[Theta][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,S,L2S,L1L2S,dSd\[Theta],dL2Sd\[Theta],dL1L2Sd\[Theta],a];
    {dfnn0d\[Omega],dfnmb0d\[Omega],dfnmb1d\[Omega],dfmbmb0d\[Omega],dfmbmb1d\[Omega],dfmbmb2d\[Omega]} = dfabid\[Omega][\[Zeta],\[Zeta]bar,sin\[Theta]p,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega],a];
    un = -((rp^2+a^2)*En - a*Lz)/(2*\[CapitalSigma]); (* Four-velocity in Kinnersley tetrad *)
    ul = -((rp^2+a^2)*En - a*Lz)/(\[CapitalDelta]);
    um = (I*(a*sin\[Theta]p^2*En - Lz) + Uz)/(-Sqrt[2]*sin\[Theta]p*\[Zeta]bar);
    umb = Conjugate[um];
    (* Dipole term *)
    Sln  = (-((rp (Kc-a^2 zp^2))/(Sqrt[Kc] \[CapitalSigma]))); (* Spin tensor in Kinnersley tetrad *)
    Snm  = (\[Zeta]/Sqrt[Kc])*um*un;
    Snmb = Conjugate[Snm];
    Slmb = (-(\[Zeta]/Sqrt[Kc]))*ul*umb;
    Smmb = ((I a zp (Kc+rp^2))/(Sqrt[Kc] \[CapitalSigma]));
    Amnn   = un^2;
    Amnmb  = un*umb;
    Ammbmb = umb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Adnn  = (Scd\[Gamma]ndc*un-Sln*2*Re[gamma]*un-2*Re[Snmb*((Conjugate[alpha]+beta)*un-mu*um)]);
    Admbmb= (Scd\[Gamma]mbdc*umb-Snmb*(-pi*ul)-Slmb*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)+Smmb*(-(-alpha+Conjugate[beta])*umb));
    Adnmb = (Scd\[Gamma]ndc*umb+Scd\[Gamma]mbdc*un-Sln*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)-Snmb*(Conjugate[rho]*un-mu*ul-(Conjugate[alpha]-beta)*umb)
      -Snm*(-(-alpha+Conjugate[beta])*umb)-Snmb*(-Conjugate[pi]*umb-pi*um)-Slmb*(2*Re[gamma]*un)+Smmb*((alpha+Conjugate[beta])*un-Conjugate[mu]*umb))/2;
    St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
    S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
    (* Derivatives of the trajectory *)
    rp1 = {correctionp[[i\[Theta]]]["\[Delta]r"], correctionp[[-i\[Theta]]]["\[Delta]r"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
    zp1 = {correctionp[[i\[Theta]]]["\[Delta]z"],-correctionp[[-i\[Theta]]]["\[Delta]z"]};
    Urp1 = {correctionp[[i\[Theta]]]["\[Delta]Ur"], correctionp[[-i\[Theta]]]["\[Delta]Ur"]};
    Uzp1 = {correctionp[[i\[Theta]]]["\[Delta]Uz"],-correctionp[[-i\[Theta]]]["\[Delta]Uz"]};
    dzpdr = {derivativesp[[i\[Theta]]]["dzdr"],-derivativesp[[-i\[Theta]]]["dzdr"]}; (* Derivatives of the coordinates and four-velocity for each quadrant *)
    dzpdx = {derivativesp[[i\[Theta]]]["dzdx"],-derivativesp[[-i\[Theta]]]["dzdx"]};
    d\[Theta]pdr = -dzpdr/sin\[Theta]p;
    d\[Theta]pdx = -dzpdx/sin\[Theta]p;
    dUzpdr = {derivativesp[[i\[Theta]]]["dUzdr"],-derivativesp[[-i\[Theta]]]["dUzdr"]};
    dUzpdx = {derivativesp[[i\[Theta]]]["dUzdx"],-derivativesp[[-i\[Theta]]]["dUzdx"]};
    \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction and derivatives of \[CapitalSigma] *)
    d\[CapitalSigma]dr = 2*(rp + a^2*zp*dzpdr);
    d\[CapitalSigma]dx = 2*a^2*zp*dzpdx;
    exp1 = I*{(\[Omega]*correctionp[[ i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[ i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
              (\[Omega]*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]t"]-m*correctionp[[-i\[Theta]]]["\[Delta]\[CapitalDelta]\[Phi]"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]}; (* Linear corrections and derivatives of the exponantial term *)
    dexpdr = I*{(\[Omega]*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]tdr"]-m*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]\[Phi]dr"])+d\[Omega]dr*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
                (\[Omega]*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]tdr"]-m*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]\[Phi]dr"])-d\[Omega]dr*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};
    dexpdx = I*{(\[Omega]*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]tdx"]-m*derivativesp[[ i\[Theta]]]["d\[CapitalDelta]\[Phi]dx"])+d\[Omega]dx*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
                (\[Omega]*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]tdx"]-m*derivativesp[[-i\[Theta]]]["d\[CapitalDelta]\[Phi]dx"])-d\[Omega]dx*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};
    un1  = -( ((2*rp*rp1)*En) - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*\[CapitalSigma]1 + 
             ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]) + Urp1))/(2*\[CapitalSigma]);(* Linear parts and derivatives of the four-velocity in Kinnersley tetrad *)
    dundr  = -( ((2*rp)*En) - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*d\[CapitalSigma]dr + 
             ((rp^2+a^2)*dEndr - a*dLzdr))/(2*\[CapitalSigma]);
    dundx  = -( - ((rp^2+a^2)*En - a*Lz)/(\[CapitalSigma])*d\[CapitalSigma]dx + 
             ((rp^2+a^2)*dEndx - a*dLzdx))/(2*\[CapitalSigma]);
    umb1 = ( (-I*(-2*a*zp*zp1*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp])) + Uzp1))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    dumbdr = ( (-I*(-2*a*zp*dzpdr*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*dzpdr/sin\[Theta]p^2 + (1-I*a*dzpdr)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*dEndr - dLzdr) + dUzpdr))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    dumbdx = ( (-I*(-2*a*zp*dzpdx*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*dzpdx/sin\[Theta]p^2 + (-I*a*dzpdx)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*dEndx - dLzdx) + dUzpdx))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*un + 2*un1)*un + Adnn; (* Source term *)
    Annt\[Phi]S = (St\[Phi]n + exp1*un)*un;
    AnnrS  = (Srn + rp1*un)*un;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*un)*un;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*un*umb + un1*umb + un*umb1 + Adnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*umb + St\[Phi]mb*un)/2 + exp1*un*umb);
    AnmbrS  = ((Srn*umb + Srmb*un)/2 + rp1*un*umb);
    Anmb\[Theta]S  = ((S\[Theta]n*umb + S\[Theta]mb*un)/2 - zp1/sin\[Theta]p*un*umb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*umb + 2*umb1)*umb + Admbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*umb)*umb;
    AmbmbrS  = (Srmb + rp1*umb)*umb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*umb)*umb;
      {Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RIn - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRIndr + Ammbmb*fmbmb2*d2RIndr2)*exp\[Theta]^{1,-1})], (* Total of all quadrants *)
       Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*(((Amnn*fnn0+Amnmb*fnmb0+Ammbmb*fmbmb0)*RUp - (Amnmb*fnmb1+Ammbmb*fmbmb1)*dRUpdr + Ammbmb*fmbmb2*d2RUpdr2)*exp\[Theta]^{1,-1})], 
       Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RIn*fnn0 + AnnrS*(dRIndr*fnn0 + RIn*dfnn0dr) + Ann\[Theta]S*RIn*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RIn*fnmb0 + AnmbrS*( RIn*dfnmb0dr +  dRIndr*fnmb0) + Anmb\[Theta]S* RIn*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRIndr*fnmb1 + AnmbrS*(dRIndr*dfnmb1dr + d2RIndr2*fnmb1) + Anmb\[Theta]S*dRIndr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RIn*fmbmb0 + AmbmbrS*(  dRIndr*fmbmb0 +   RIn*dfmbmb0dr) + Ambmb\[Theta]S*  RIn*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRIndr*fmbmb1 + AmbmbrS*( d2RIndr2*fmbmb1 +  dRIndr*dfmbmb1dr) + Ambmb\[Theta]S* dRIndr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RIndr2*fmbmb2 + AmbmbrS*(d3RIndr3*fmbmb2 + d2RIndr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RIndr2*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        dw\[Theta]du\[Theta]*\[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*dRInd\[Omega] - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*d2RIndrd\[Omega] + Ammbmb*fmbmb2*d3RIndr2d\[Omega]
          + (Amnn*dfnn0d\[Omega] + Amnmb*dfnmb0d\[Omega] + Ammbmb*dfmbmb0d\[Omega])*RIn - (Amnmb*dfnmb1d\[Omega] + Ammbmb*dfmbmb1d\[Omega])*dRIndr + Ammbmb*dfmbmb2d\[Omega]*d2RIndr2)*\[Omega]1*exp\[Theta]^{1,-1})], 
       Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUp*fnn0 + AnnrS*(dRUpdr*fnn0 + RUp*dfnn0dr) + Ann\[Theta]S*RUp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUp*fnmb0 + AnmbrS*( RUp*dfnmb0dr +  dRUpdr*fnmb0) + Anmb\[Theta]S* RUp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUpdr*fnmb1 + AnmbrS*(dRUpdr*dfnmb1dr + d2RUpdr2*fnmb1) + Anmb\[Theta]S*dRUpdr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUp*fmbmb0 + AmbmbrS*(  dRUpdr*fmbmb0 +   RUp*dfmbmb0dr) + Ambmb\[Theta]S*  RUp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUpdr*fmbmb1 + AmbmbrS*( d2RUpdr2*fmbmb1 +  dRUpdr*dfmbmb1dr) + Ambmb\[Theta]S* dRUpdr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RUpdr2*fmbmb2 + AmbmbrS*(d3RUpdr3*fmbmb2 + d2RUpdr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RUpdr2*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        dw\[Theta]du\[Theta]*\[CapitalSigma]*(((Amnn*fnn0 + Amnmb*fnmb0 + Ammbmb*fmbmb0)*dRUpd\[Omega] - (Amnmb*fnmb1 + Ammbmb*fmbmb1)*d2RUpdrd\[Omega] + Ammbmb*fmbmb2*d3RUpdr2d\[Omega]
          + (Amnn*dfnn0d\[Omega] + Amnmb*dfnmb0d\[Omega] + Ammbmb*dfmbmb0d\[Omega])*RUp - (Amnmb*dfnmb1d\[Omega] + Ammbmb*dfmbmb1d\[Omega])*dRUpdr + Ammbmb*dfmbmb2d\[Omega]*d2RUpdr2)*\[Omega]1*exp\[Theta]^{1,-1})], 
      Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Amnn*fnn0*RIn + Amnmb*(fnmb0*RIn - fnmb1*dRIndr) + Ammbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundr*fnn0*RIn + (un*dumbdr+dundr*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdr*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Amnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RIn + fnn0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)) + 
               Amnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RIn + fnmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRIndr - fnmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr)) + 
              Ammbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RIn + fmbmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRIndr - fmbmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RIndr2 + fmbmb2*(d3RIndr3 + d3RIndr2d\[Omega]*d\[Omega]dr)))) 
              )*exp\[Theta]^{1,-1}], 
       Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Amnn*fnn0*RUp + Amnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Ammbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundr*fnn0*RUp + (un*dumbdr+dundr*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdr*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Amnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RUp + fnn0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)) + 
               Amnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RUp + fnmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fnmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr)) + 
              Ammbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RUp + fmbmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fmbmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RUpdr2 + fmbmb2*(d3RUpdr3 + d3RUpdr2d\[Omega]*d\[Omega]dr)))) 
              )*exp\[Theta]^{1,-1}], 
      Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Amnn*fnn0*RIn + Amnmb*(fnmb0*RIn - fnmb1*dRIndr) + Ammbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundx*fnn0*RIn + (un*dumbdx+dundx*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdx*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Amnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RIn + fnn0*(dRInd\[Omega]*d\[Omega]dx)) + 
               Amnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RIn + fnmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRIndr - fnmb1*(d2RIndrd\[Omega]*d\[Omega]dx)) + 
              Ammbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RIn + fmbmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfmbmb1d\[Theta]*d\[Theta]pdx + dfmbmb1d\[Omega]*d\[Omega]dx)*dRIndr - fmbmb1*(d2RIndrd\[Omega]*d\[Omega]dx) + 
                      (dfmbmb2d\[Theta]*d\[Theta]pdx + dfmbmb2d\[Omega]*d\[Omega]dx)*d2RIndr2 + fmbmb2*(d3RIndr2d\[Omega]*d\[Omega]dx)))) 
              )*exp\[Theta]^{1,-1}], 
       Total[dw\[Theta]du\[Theta]*\[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Amnn*fnn0*RUp + Amnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Ammbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundx*fnn0*RUp + (un*dumbdx+dundx*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdx*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Amnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RUp + fnn0*(dRUpd\[Omega]*d\[Omega]dx)) + 
               Amnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RUp + fnmb0*(dRUpd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRUpdr - fnmb1*(d2RUpdrd\[Omega]*d\[Omega]dx)) + 
              Ammbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RUp + fmbmb0*(dRUpd\[Omega]*d\[Omega]dx)
                    - (dfmbmb1d\[Theta]*d\[Theta]pdx + dfmbmb1d\[Omega]*d\[Omega]dx)*dRUpdr - fmbmb1*(d2RUpdrd\[Omega]*d\[Omega]dx) + 
                      (dfmbmb2d\[Theta]*d\[Theta]pdx + dfmbmb2d\[Omega]*d\[Omega]dx)*d2RUpdr2 + fmbmb2*(d3RUpdr2d\[Omega]*d\[Omega]dx)))) 
              )*exp\[Theta]^{1,-1}]
       },
    {i\[Theta], 1, steps\[Theta]/2}
  ];
  W = (RIn*dRUpdr - dRIndr*RUp)/\[CapitalDelta]; (* Invariant Wronskian *)
  dWd\[Omega]= (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/\[CapitalDelta]; (* \[Omega]-derivative of the invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalGamma]*W*steps\[Theta]); (* Geodesic amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalGamma]*W*steps\[Theta]);
  CPlus1  = 2*Pi*(sumPlus1  - \[CapitalGamma]1/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*\[Omega]1/W*sumPlus0 )/(\[CapitalGamma]*W*steps\[Theta]); (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1 - \[CapitalGamma]1/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*\[Omega]1/W*sumMinus0)/(\[CapitalGamma]*W*steps\[Theta]);
  dCPlusdr  = 2*Pi*(dsumPlusdr  - d\[CapitalGamma]dr/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*d\[Omega]dr/W*sumPlus0 )/(\[CapitalGamma]*W*steps\[Theta]); (* Derivatives of the amplitudes *)
  dCMinusdr = 2*Pi*(dsumMinusdr - d\[CapitalGamma]dr/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*d\[Omega]dr/W*sumMinus0)/(\[CapitalGamma]*W*steps\[Theta]);
  dCPlusdx  = 2*Pi*(dsumPlusdx  - d\[CapitalGamma]dx/\[CapitalGamma]*sumPlus0  - dWd\[Omega]*d\[Omega]dx/W*sumPlus0 )/(\[CapitalGamma]*W*steps\[Theta]);
  dCMinusdx = 2*Pi*(dsumMinusdx - d\[CapitalGamma]dx/\[CapitalGamma]*sumMinus0 - dWd\[Omega]*d\[Omega]dx/W*sumMinus0)/(\[CapitalGamma]*W*steps\[Theta]);
  (*Print[NumberForm[Re[CPlus0],16]];*)
  <|
    "l" -> l,
    "m" -> m,
    "k" -> k,
    "n" -> 0,
    "\[Omega]" -> \[Omega],
    "\[Omega]Correction" -> \[Omega]1,
    "\[Omega]Derivatives" -> <|
      "r" -> d\[Omega]dr,
      "x" -> d\[Omega]dx
    |>,
    "Amplitudes" -> <|
      "\[ScriptCapitalI]" -> CPlus0,
      "\[ScriptCapitalH]" -> CMinus0
    |>,
    "AmplitudesCorrection" -> <|
      "\[ScriptCapitalI]" -> CPlus1,
      "\[ScriptCapitalH]" -> CMinus1
    |>,
    "AmplitudesDerivatives" -> <|
      "r" -> <|
        "\[ScriptCapitalI]" -> dCPlusdr,
        "\[ScriptCapitalH]" -> dCMinusdr
      |>,
      "x" -> <|
        "\[ScriptCapitalI]" -> dCPlusdx,
        "\[ScriptCapitalH]" -> dCMinusdx
      |>
    |>,
    "\[Alpha]" -> \[Alpha],
    "d\[Alpha]d\[Omega]" -> d\[Alpha]d\[Omega],
    "S" -> N[SWSH[Pi/2,0][[1]]],
    "Fluxes" -> <|
      "Energy" -> <|
        "\[ScriptCapitalI]" -> Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]" -> \[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum" -> <|
        "\[ScriptCapitalI]" -> Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]" -> \[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection" -> <|
      "Energy" -> <|
        "\[ScriptCapitalI]" -> (2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]" -> \[Alpha]*(2*Re[CMinus1*Conjugate[CMinus0]] + (d\[Alpha]d\[Omega]/\[Alpha]*Abs[CMinus0]^2 - 2*Abs[CMinus0]^2/\[Omega])*\[Omega]1)/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum" -> <|
        "\[ScriptCapitalI]" -> (2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>
    |>,
    "FluxesDerivatives" -> <|
      "r" -> <|
        "Energy"  ->  <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdr*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*d\[Omega]dr/\[Omega])/(4Pi*\[Omega]^2),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dr/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdr*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*d\[Omega]dr/\[Omega])/(4Pi*\[Omega]^2)
        |>,
        "AngularMomentum" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdr*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*d\[Omega]dr/\[Omega])*m/(4Pi*\[Omega]^3),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dr/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdr*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*d\[Omega]dr/\[Omega])*m/(4Pi*\[Omega]^3)
        |>
      |>,
      "x" -> <|
        "Energy" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdx*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*d\[Omega]dx/\[Omega])/(4Pi*\[Omega]^2),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dx/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdx*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*d\[Omega]dx/\[Omega])/(4Pi*\[Omega]^2)
        |>,
        "AngularMomentum" -> <|
          "\[ScriptCapitalI]" -> (2*Re[dCPlusdx*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*d\[Omega]dx/\[Omega])*m/(4Pi*\[Omega]^3),
          "\[ScriptCapitalH]" -> \[Alpha]*(d\[Alpha]d\[Omega]*d\[Omega]dx/\[Alpha]*Abs[CMinus0]^2 + 2*Re[dCMinusdx*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*d\[Omega]dx/\[Omega])*m/(4Pi*\[Omega]^3)
        |>
      |>
    |>,
    "\[Lambda]" -> \[Lambda],
    "d\[Lambda]d\[Omega]" -> d\[Lambda]d\[Omega],
    "steps\[Theta]" -> steps\[Theta]
  |>
]


Options[TeukolskySpinModeSphericalCorrectionAnalytical] = {WorkingPrecision->32};
TeukolskySpinModeSphericalCorrectionAnalytical[l_?IntegerQ,m_?IntegerQ,k_?IntegerQ,orbit_,{angparNew_,TeukolskySolverHS1spin_},OptionsPattern[]]:=Module[{
    a,p,e,x,En0,Lz0,K0,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi],z,\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t1,\[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]z,\[CapitalUpsilon]\[Tau]z,\[Omega],\[Omega]prec,\[Omega]1,
    SWSH,dSWSHd\[Omega],R,\[Lambda],d\[Lambda]d\[Omega],\[ScriptCapitalC]2,d\[ScriptCapitalC]2d\[Omega],rplus,P,\[Epsilon],\[Alpha],d\[Alpha]d\[Omega],
    sumPlus0,sumPlus1,sumMinus0,sumMinus1,stepsz,iz,wz,
    rp,\[CapitalDelta],d\[CapitalDelta],K,dKdr,dKd\[Omega],d2Kdrd\[Omega],V,dVd\[Omega],
    Rrp,RIn,dRIndr,d2RIndr2,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    DRIn,dDRIndr,DDRIn,DRUp,dDRUpdr,DDRUp,dDRInd\[Omega],d2DRIndrd\[Omega],dDDRInd\[Omega],dDRUpd\[Omega],d2DRUpdrd\[Omega],dDDRUpd\[Omega],
    zp,\[Theta]p,Uzp,expz,expz1,sin\[Theta]p,K\[Theta],dK\[Theta]d\[Omega],S,dSd\[Theta],d2Sd\[Theta]2,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],
    L2S,dL2Sd\[Theta],L1L2S,dL2Sd\[Omega],d2L2Sd\[Theta]d\[Omega],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,vn,vmb,
    FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn,dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega],
    FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp,dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega],
    W,W1,CPlus0,CPlus1,CMinus0,CMinus1},
  If[l < 2 || Abs[m] > l, Return[$Failed]];
  a = orbit["a"];(* Orbital parameters *)
  p = orbit["p"];
  e = orbit["e"];
  x = orbit["Inclination"];
  En0 = orbit["Energy"]; (* Shifted constants of motion *)
  Lz0 = orbit["AngularMomentum"];
  K0 = orbit["CarterConstant"] + (Lz0 - a*En0)^2;
  {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[orbit["Frequencies"]];
  {\[CapitalOmega]r,\[CapitalOmega]z,\[CapitalOmega]\[Phi]} = {\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t; (* BL frequencies *)
  z[qz_] := Cos[orbit["Trajectory"][[3]][qz]];
  \[CapitalUpsilon]\[Tau] = orbit["ProperTimeFrequency"];
  \[CapitalUpsilon]t1 = -3*\[CapitalUpsilon]\[Tau]/(2*Sqrt[K0]);
  \[CapitalUpsilon]tz = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencyt\[Theta][a,p,e,x,{En0,Lz0,K0 - (Lz0 - a*En0)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,p,e,x]] + a*Lz0;
  \[CapitalUpsilon]\[Phi]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Phi]\[Theta][a,p,e,x,{En0,Lz0,K0 - (Lz0 - a*En0)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,p,e,x]] - a*En0;
  \[CapitalUpsilon]\[Tau]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Tau]z[a,p,e,x];
  Print["\[CapitalUpsilon]\[Tau]z = "<>ToString[\[CapitalUpsilon]\[Tau]z]];
  \[Omega] = m*\[CapitalOmega]\[Phi] + k*\[CapitalOmega]z; (* Frequency of mode *)
  \[Omega]prec = m*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"]
        + k*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)"]; (* Frequency with higher precision *)
  \[Omega]1 = 3*\[CapitalUpsilon]\[Tau]*\[Omega]/(2*Sqrt[K0]*\[CapitalUpsilon]t);
  If[!(\[Omega]\[Element]Reals), Return[$Failed]];
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]} = angparNew[-2,l,m,
                                     SetPrecision[a,    OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],
                                     SetPrecision[\[Omega]prec,OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)],1,
                                           "precODE" -> OptionValue[WorkingPrecision]+(5+1.5*\[Omega]+5.3*(l-Max[2,Abs[m]])+2.6*a)];(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  (*Print["Calculating radial functions"<>ToString@AbsoluteTiming[*)R = TeukolskySolverHS1spin[p,-2,l,m,
                                                                SetPrecision[a,    OptionValue[WorkingPrecision]+5],
                                                                SetPrecision[\[Omega]prec,OptionValue[WorkingPrecision]+5],1,
                                                                SetPrecision[\[Lambda],    OptionValue[WorkingPrecision]+5],
                                                                SetPrecision[d\[Lambda]d\[Omega], OptionValue[WorkingPrecision]+5],
                                                                       "precODE" -> OptionValue[WorkingPrecision]](*][[1]]]*);
  \[ScriptCapitalC]2 = ((\[Lambda]+2)^2+4a*\[Omega](m-a*\[Omega]))*(\[Lambda]^2+36a*\[Omega](m-a*\[Omega]))-(2\[Lambda]+3)*(48a*\[Omega](m-2a*\[Omega]))+144*\[Omega]^2*(1-a^2); (*  TS constant *)
  d\[ScriptCapitalC]2d\[Omega] = 4 \[Lambda]^3 d\[Lambda]d\[Omega]+4 \[Lambda]^2 (3 d\[Lambda]d\[Omega] + 10 a (m-2 a \[Omega]))+8 \[Lambda] (d\[Lambda]d\[Omega] (1+10 a m \[Omega]-10 a^2 \[Omega]^2)+6 a (m+2 a \[Omega])) + 
        48 \[Omega] (a m d\[Lambda]d\[Omega]+6-18 a^3 m \[Omega]+12 a^4 \[Omega]^2+a^2 (d\[Lambda]d\[Omega] \[Omega]+6 m^2));  (* \[Omega] derivative of the TS constant *)
  rplus = 1+Sqrt[1-a^2];  (*  horizon r_+  *)
  P = \[Omega]-m*a/(2*rplus); (* frequency at the horizon *)
  \[Epsilon] = Sqrt[1^2-a^2]/(4*rplus);
  \[Alpha] = 256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2; (* constant for horizon fluxes *)
  d\[Alpha]d\[Omega] = -256*(2*rplus)^5*P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3/\[ScriptCapitalC]2^2*d\[ScriptCapitalC]2d\[Omega] + 256*(2*rplus)^5*((P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 +
       P*(2*P)*(P^2+16*\[Epsilon]^2)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(2*P)*\[Omega]^3 + P*(P^2+4*\[Epsilon]^2)*(P^2+16*\[Epsilon]^2)*3*\[Omega]^2)/\[ScriptCapitalC]2; (* \[Omega] derivative of the constant for horizon fluxes *)
  (* numbers of steps for wr and wz integration *)
  stepsz = Max[32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               32*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsz]<>" steps in wz"];
  rp = p;
  \[CapitalDelta]  = rp^2-2rp+a^2;
  K  = (rp^2+a^2)*\[Omega]-a*m;
  dKdr = 2*rp*\[Omega];
  dKd\[Omega]  = (rp^2+a^2);
  d2Kdrd\[Omega] = 2*rp;
  d\[CapitalDelta] = 2*(rp-1);
  V  = -(K^2 + 4I*(rp-1)*K)/\[CapitalDelta] + 8*I*\[Omega]*rp + \[Lambda]; (* Potential in radial Teukolsky equation *)
  dVd\[Omega]  = -(2*K*dKd\[Omega] + 4I*(rp-1)*dKd\[Omega])/\[CapitalDelta] + 8*I*rp + d\[Lambda]d\[Omega]; (* \[Omega]-derivative of the potential in radial Teukolsky equation *)
  {{{RIn,dRIndr},{dRInd\[Omega],d2RIndrd\[Omega]}},{{RUp,dRUpdr},{dRUpd\[Omega],d2RUpdrd\[Omega]}}} = R; (* Radial function and its derivatives *)
  d2RIndr2  = (V*RIn + d\[CapitalDelta]*dRIndr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  d3RIndr2d\[Omega] = (V*dRInd\[Omega] + dVd\[Omega]*RIn + d\[CapitalDelta]*d2RIndrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRIn = dRIndr - I*K/\[CapitalDelta]*RIn;
  dDRInd\[Omega] = d2RIndrd\[Omega] - I*K/\[CapitalDelta]*dRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RIn;
  dDRIndr = d2RIndr2 - I*K/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  d2DRIndrd\[Omega] = d3RIndr2d\[Omega] - I*K/\[CapitalDelta]*d2RIndrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRIndr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRInd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RIn;
  DDRIn = dDRIndr - I*K/\[CapitalDelta]*DRIn;
  dDDRInd\[Omega] = d2DRIndrd\[Omega] - I*K/\[CapitalDelta]*dDRInd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRIn;
  (* Up solutions of the radial equation and their r and \[Omega] derivatives *)
  d2RUpdr2  = (V*RUp + d\[CapitalDelta]*dRUpdr)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  d3RUpdr2d\[Omega] = (V*dRUpd\[Omega] + dVd\[Omega]*RUp + d\[CapitalDelta]*d2RUpdrd\[Omega])/\[CapitalDelta];
  (* Operators D and their r and \[Omega] derivatives *)
  DRUp = dRUpdr - I*K/\[CapitalDelta]*RUp;
  dDRUpd\[Omega] = d2RUpdrd\[Omega] - I*K/\[CapitalDelta]*dRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*RUp;
  dDRUpdr = d2RUpdr2 - I*K/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  d2DRUpdrd\[Omega] = d3RUpdr2d\[Omega] - I*K/\[CapitalDelta]*d2RUpdrd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*dRUpdr - I*( dKdr*\[CapitalDelta] - K*d\[CapitalDelta])/\[CapitalDelta]^2*dRUpd\[Omega] - I*( d2Kdrd\[Omega]*\[CapitalDelta] - dKd\[Omega]*d\[CapitalDelta])/\[CapitalDelta]^2*RUp;
  DDRUp = dDRUpdr - I*K/\[CapitalDelta]*DRUp;
  dDDRUpd\[Omega] = d2DRUpdrd\[Omega] - I*K/\[CapitalDelta]*dDRUpd\[Omega] - I*dKd\[Omega]/\[CapitalDelta]*DRUp;
  (* integration over w\[Theta] *)
  {sumPlus0, sumPlus1, sumMinus0, sumMinus1} = Sum[
    wz = N[(iz-1/2)*2Pi/stepsz,Precision[{a,p,e,x}]];
    zp = z[wz];
    \[Theta]p = ArcCos[zp];
    Uzp = {1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(K0-a^2*zp^2)];(* Polar geodesic velocity *)
    expz = Exp[I*(\[Omega]*(orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][wz])-m*(orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][wz])+2Pi*k*(iz-1/2)/stepsz)];(* Exponential term with geodesic \[CapitalDelta]tz and \[CapitalDelta]\[Phi]z *)
    expz1 = {1,-1}*3*I*\[Omega]/(2*Sqrt[K0])*(\[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][wz] - orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Tau]\[Theta]"][wz]);(* Linear part of the exponential term *)
    sin\[Theta]p = Sqrt[1-zp^2];
    K\[Theta] = a*\[Omega]*(1-zp^2) - m;
    dK\[Theta]d\[Omega] = a*(1-zp^2);
    (* Spin-weighted spheroidal harmonics and their \[Theta]-derivatives *)
    {S,dSd\[Theta]}=SWSH[\[Theta]p]; 
    d2Sd\[Theta]2 = -(zp/sin\[Theta]p)*dSd\[Theta] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*S;
    (* \[Omega]-derivatives of S *)
    {dSd\[Omega],d2Sd\[Theta]d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
    d3Sd\[Theta]2d\[Omega] = -(zp/sin\[Theta]p)*d2Sd\[Theta]d\[Omega] - (-a^2*\[Omega]^2*(1-zp^2) - (m-2*zp)^2/(1-zp^2) + 4*a*\[Omega]*zp - 2 + 2*m*a*\[Omega] + \[Lambda])*dSd\[Omega] - (-2*a^2*\[Omega]*(1-zp^2) + 4*a*zp + 2*m*a + d\[Lambda]d\[Omega])*S;
    (* Operators L and their \[Omega]-derivatives *)
    L2S = dSd\[Theta] + (K\[Theta] + 2 zp)*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
    dL2Sd\[Theta] = d2Sd\[Theta]2 + (K\[Theta] + 2 zp)*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*S;
    L1L2S = dL2Sd\[Theta] + (K\[Theta] + zp)*L2S/sin\[Theta]p; 
    dL2Sd\[Omega] = d2Sd\[Theta]d\[Omega] + (K\[Theta] + 2 zp)*dSd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*S/sin\[Theta]p;(* Operators acting on S(\[Theta]) and derivatives of these operators *)
    d2L2Sd\[Theta]d\[Omega] = d3Sd\[Theta]2d\[Omega] + (K\[Theta] + 2 zp)*d2Sd\[Theta]d\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*dSd\[Theta]/sin\[Theta]p + (2*(a*\[Omega]*zp - 1) - (K\[Theta] + 2 zp)/sin\[Theta]p^2*zp)*dSd\[Omega] + (2*a*zp - dK\[Theta]d\[Omega]/sin\[Theta]p^2*zp)*S;
    dL1L2Sd\[Omega] = d2L2Sd\[Theta]d\[Omega] + (K\[Theta] + zp)*dL2Sd\[Omega]/sin\[Theta]p + dK\[Theta]d\[Omega]*L2S/sin\[Theta]p; 
    (* Integration f functions of r and \[Theta] *)
    \[Zeta] = rp-I*a*zp;
    \[Zeta]bar = rp+I*a*zp;
    vn = -((rp^2+a^2)*En0 - a*Lz0 )/\[CapitalDelta]; (* Four-velocity in rotated Kinnersley tetrad *)
    vmb = (-I*(a*sin\[Theta]p^2*En0 - Lz0) + Uzp)/sin\[Theta]p^2;
    (* In solution *)
    {FnnIn,FnmbIn,FmbmbIn,GnIn,GmbIn} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RIn,DRIn,DDRIn,S,L2S,L1L2S];(* F_ab functions *)
    {dFnnInd\[Omega],dFnmbInd\[Omega],dFmbmbInd\[Omega]}   = dFabd\[Omega][\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RIn,DRIn,DDRIn,dRInd\[Omega],dDRInd\[Omega],dDDRInd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
    (* Up solution *)
    {FnnUp,FnmbUp,FmbmbUp,GnUp,GmbUp} = Fab[\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RUp,DRUp,DDRUp,S,L2S,L1L2S];(* F_ab functions *)
    {dFnnUpd\[Omega],dFnmbUpd\[Omega],dFmbmbUpd\[Omega]}   = dFabd\[Omega][\[Zeta],\[Zeta]bar,a,sin\[Theta]p,RUp,DRUp,DDRUp,dRUpd\[Omega],dDRUpd\[Omega],dDDRUpd\[Omega],S,L2S,L1L2S,dSd\[Omega],dL2Sd\[Omega],dL1L2Sd\[Omega]];(* \[Omega]-derivatives of F_ab *)
    (* Summing the geodesic ampltude over all quadrants *)
    (* Summing the geodesic ampltude over all quadrants *)
    {
      Total[(vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*expz^{1,-1}],
      Total[(
        (vn*( GnIn ) + vmb*( GmbIn ))/Sqrt[K0] +
        (vn*vn*FnnIn + vn*vmb*FnmbIn + vmb*vmb*FmbmbIn)*(expz1) +
        (vn*vn*dFnnInd\[Omega] + vn*vmb*dFnmbInd\[Omega] + vmb*vmb*dFmbmbInd\[Omega])*\[Omega]1
       )*expz^{1,-1}],
      Total[(vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*expz^{1,-1}],
      Total[(
        (vn*( GnUp ) + vmb*( GmbUp ))/Sqrt[K0] +
        (vn*vn*FnnUp + vn*vmb*FnmbUp + vmb*vmb*FmbmbUp)*(expz1) +
        (vn*vn*dFnnUpd\[Omega] + vn*vmb*dFnmbUpd\[Omega] + vmb*vmb*dFmbmbUpd\[Omega])*\[Omega]1
      )*expz^{1,-1}]
    },
    {iz, 1, stepsz/2}];
  W = (RIn*dRUpdr - dRIndr*RUp)/(rp^2-2*rp+a^2); (* Invariant Wronskian *)
  W1 = (dRInd\[Omega]*dRUpdr + RIn*d2RUpdrd\[Omega] - d2RIndrd\[Omega]*RUp - dRIndr*dRUpd\[Omega])/(rp^2-2*rp+a^2)*\[Omega]1; (* Invariant Wronskian *)
  CPlus0  = 2*Pi*sumPlus0/(\[CapitalUpsilon]t*W*stepsz); (* Amplitudes *)
  CMinus0 = 2*Pi*sumMinus0/(\[CapitalUpsilon]t*W*stepsz);
  CPlus1  = 2*Pi*(sumPlus1)/(\[CapitalUpsilon]t*W*stepsz) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CPlus0; (* Linear parts of the amplitudes *)
  CMinus1 = 2*Pi*(sumMinus1)/(\[CapitalUpsilon]t*W*stepsz) + (En0/(2*Sqrt[K0]) - \[CapitalUpsilon]t1/\[CapitalUpsilon]t - W1/W)*CMinus0;
  <|
    "l"->l,
    "m"->m,
    "k"->k,
    "n"->0,
    "\[Omega]"->\[Omega],
    "\[Omega]Correction"->\[Omega]1,
    "Amplitudes"->
    <|
      "\[ScriptCapitalI]"->CPlus0,
      "\[ScriptCapitalH]"->CMinus0
    |>,
    "AmplitudesCorrection"->
    <|
      "\[ScriptCapitalI]"->CPlus1,
      "\[ScriptCapitalH]"->CMinus1
    |>,
    "\[Alpha]"->\[Alpha],
    "S"->SWSH[Pi/2,0],
    "Fluxes"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->Abs[CPlus0]^2*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*Abs[CMinus0]^2*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"->(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*Abs[CPlus0]^2/(2Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*\[Alpha]*Abs[CMinus0]^2/(2Pi*\[Omega]^3)
      |>
    |>,
    "FluxesCorrection"->
    <|
      "Energy"-><|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 2*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2),
        "\[ScriptCapitalH]"->\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 2*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(4Pi*\[Omega]^2)
      |>,
      "AngularMomentum"->
      <|
        "\[ScriptCapitalI]"->(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3),
        "\[ScriptCapitalH]"->\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])*m/(4Pi*\[Omega]^3)
      |>,
      "CarterConstantK"->
      <|
        "\[ScriptCapitalI]"-> (En0*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z) - a*( m - a*\[Omega] ) - 3*\[Omega]*( -\[CapitalUpsilon]\[Tau]z + \[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*\[CapitalUpsilon]tz ))/(2*Sqrt[K0])*Abs[CPlus0]^2/(2Pi*\[Omega]^3)+(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*(2*Re[CPlus1*Conjugate[CPlus0]] - 3*Abs[CPlus0]^2*\[Omega]1/\[Omega])/(2Pi*\[Omega]^3),
        "\[ScriptCapitalH]"-> (En0*(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z) - a*( m - a*\[Omega] ) - 3*\[Omega]*( -\[CapitalUpsilon]\[Tau]z + \[CapitalUpsilon]\[Tau]/\[CapitalUpsilon]t*\[CapitalUpsilon]tz ))/(2*Sqrt[K0])*\[Alpha]*Abs[CMinus0]^2/(2Pi*\[Omega]^3)+(k*\[CapitalUpsilon]z - \[Omega]*\[CapitalUpsilon]tz + m*\[CapitalUpsilon]\[Phi]z)*\[Alpha]*(d\[Alpha]d\[Omega]*\[Omega]1/\[Alpha]*Abs[CMinus0]^2 + 2*Re[CMinus1*Conjugate[CMinus0]] - 3*Abs[CMinus0]^2*\[Omega]1/\[Omega])/(2Pi*\[Omega]^3)
      |>
    |>,
    "steps\[Theta]"->stepsz
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];

(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["TeukolskySpinFluxes`",
  {"KerrGeodesics`",
  (*"Teukolsky`",*)
  "SpinningOrbit`"}
];


TeukolskySpinMode::usage = "TeukolskySpinMode[l,m,n,k,orbit] calculates Teukolsky amplitudes and fluxes from a spinning orbit";
TeukolskySpinMode2::usage = "TeukolskySpinMode2[l,m,n,k,orbitCorrection,s] calculates Teukolsky amplitudes and fluxes from a correction to an orbit and a value of the spin";
TeukolskySpinModeCorrectionNum::usage = "TeukolskySpinModeCorrectionNum[l,m,n,k,orbitCorrection,\[Delta]\[Omega]] calculates linear correction to the fluxes from numerical derivatives of R ans S with step \[Delta]\[Omega]";
TeukolskySpinModeCorrection::usage = "TeukolskySpinModeCorrection[l,m,n,k,orbitCorrection,{angparNew,RCorrection}] calculates linear correction to the fluxes using given function for the R and S derivatives";
TeukolskySpinModeSpherical::usage = "TeukolskySpinModeSpherical[l,m,k,orbitCorrection,s] calculates Teukolsky amplitudes and fluxes from a correction to a spherical orbit and a value of the spin";
TeukolskySpinModeSphericalCorrectionNum::usage = "TeukolskySpinModeSphericalCorrectionNum[l,m,k,orbitCorrection,\[Delta]\[Omega]] calculates linear correction to the fluxes for spherical orbit from numerical derivatives of R ans S with step \[Delta]\[Omega]";
TeukolskySpinModeSphericalCorrection::usage = "TeukolskySpinModeSphericalCorrection[l,m,k,orbitCorrection,{angparNew,TeukolskySolverHS1spin}] calculates linear correction to the fluxes for spherical orbit using given function for the R and S derivatives";
TeukolskySpinModeSphericalCorrectionNew::usage = "TeukolskySpinModeSphericalCorrectionNew[l,m,k,orbitCorrection,{angparNew,TeukolskySolverHS1spin}] calculates linear correction to the fluxes for spherical orbit using given function for the R and S derivatives";
TeukolskySpinModeSphericalCorrection2::usage = "TeukolskySpinModeSphericalCorrection2[l,m,k,orbitCorrection,{angparNew,TeukolskySolverHS1spin}] calculates linear correction to the fluxes for spherical orbit using given function for the R and S derivatives";


Begin["`Private`"];


(* ::Section:: *)
(*Spinning fluxes*)


(* ::Subsection::Closed:: *)
(*Functions f_ab^(i)*)


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



(* ::Subsection::Closed:: *)
(*Circular*)


Options[TeukolskySpinModeCircularCorrection] = {WorkingPrecision->30};

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
]


(* ::Subsection::Closed:: *)
(*Generic*)


TeukolskySpinMode[l_,m_,n_,k_,orbit_]:=Module[{a,p,e,\[ScriptCapitalI],s,En0,Lz0,Kc0,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalGamma],\[Omega],SWSH,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Theta]list,
    R,\[Lambda],W,sumPlus,sumMinus,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Urp,Uzp,udtp,ud\[Phi]p,Kp,\[CapitalDelta],dK,d\[CapitalDelta],\[CapitalSigma],\[Zeta],\[Zeta]bar,fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,
    dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,C0nn,C0nmb,C0mbmb,
    Crnn,Crnmb,Crmbmb,C\[Theta]nn,C\[Theta]nmb,C\[Theta]mbmb,Cmnn,Cmnmb,Cmmbmb,Cdnn,Cdnmb,Cdmbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,
    A0,A1,A2,B1,B2,B3,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,RUprp,dRUprp,ddRUprp,dddRUprp,CPlus,CMinus,\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],expr,exp\[Theta],exp1,stepsr,steps\[Theta],h1,h2,h3,
    r0,z0,Ur0,Uz0,r0p,z0p,\[CapitalDelta]0,d\[CapitalDelta]0,K0,dK0,r1p,z1p,sin\[Theta]0p,\[Zeta]0,\[Zeta]0bar,Ur1p,Uz1p,correction,correctionp},
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Kc0 (r^2-3 a^2 z^2)))/(Sqrt[Kc0] (r^2+a^2 z^2)^3);
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
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
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
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn   = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb = (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb  = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)-Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*Kp/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
      St\[Phi]mb = -I*Kp*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
      Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
      Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
      S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
      S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
      C0nn   = (Cmnn+Cdnn+St\[Phi]n*vn);
      Crnn   = (Cmnn*s*r1p+Srn*vn);
      C\[Theta]nn   = (-Cmnn*s*z1p/sin\[Theta]0p+S\[Theta]n*vn);
      C0nmb  = (Cmnmb+Cdnmb+(St\[Phi]n*vmb+St\[Phi]mb*vn)/2);
      Crnmb  = (Cmnmb*s*r1p+(Srn*vmb+Srmb*vn)/2);
      C\[Theta]nmb  = (-Cmnmb*s*z1p/sin\[Theta]0p+(S\[Theta]n*vmb+S\[Theta]mb*vn)/2);
      C0mbmb = (Cmmbmb+Cdmbmb+St\[Phi]mb*vmb);
      Crmbmb = (Cmmbmb*s*r1p+Srmb*vmb);
      C\[Theta]mbmb = (-Cmmbmb*s*z1p/sin\[Theta]0p+S\[Theta]mb*vmb);
      A0 = C0nn*fnn0+Crnn*dfnn0dr+C\[Theta]nn*dfnn0d\[Theta]+C0nmb*fnmb0+Crnmb*dfnmb0dr+C\[Theta]nmb*dfnmb0d\[Theta]+C0mbmb*fmbmb0+Crmbmb*dfmbmb0dr+C\[Theta]mbmb*dfmbmb0d\[Theta];
      A1 = C0nmb*fnmb1+Crnmb*dfnmb1dr+C\[Theta]nmb*dfnmb1d\[Theta]+C0mbmb*fmbmb1+Crmbmb*dfmbmb1dr+C\[Theta]mbmb*dfmbmb1d\[Theta];
      A2 = C0mbmb*fmbmb2+Crmbmb*dfmbmb2dr+C\[Theta]mbmb*dfmbmb2d\[Theta];
      B1 = -Srn*vn*fnn0-(Srn*vmb+Srmb*vn)/2*fnmb0-Srmb*vmb*fmbmb0;
      B2 = -(Srn*vmb+Srmb*vn)/2*fnmb1-Srmb*vmb*fmbmb1;
      B3 = -Srmb*vmb*fmbmb2;
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


TeukolskySpinMode2[l_,m_,n_,k_,orbitCorrection_,s_]:=Module[{h1,h2,h3,a,p,e,\[ScriptCapitalI],En0,Lz0,Kc0,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],correction,r,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,
    sumPlus,sumMinus,stepsr,steps\[Theta],\[Theta]list,correctionp,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,
    RUprp,dRUprp,ddRUprp,dddRUprp,\[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,
    dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,
    Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,r1p,z1p,Ur1p,Uz1p,\[CapitalSigma]1,exp1,vn1,vmb1,
    Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus,CMinus},
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc0 (r^2-3 a^2 z^2)))/(Sqrt[Kc0] (r^2+a^2 z^2)^3);
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
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  r = orbitCorrection["TrajectoryGeo"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega] = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Frequency *)
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
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[stepsr]<>" steps in wr, "<>ToString[steps\[Theta]]<>" steps in w\[Theta]"];
  \[Theta]list = {};(* List for functions of \[Theta] *)
  correctionp = Table[correction[N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]],N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,\[ScriptCapitalI]}]]],{ir,1,stepsr/2},{i\[Theta],1,steps\[Theta]/2}];(* corrections to the trajectory at all points *)
  For[ir = 1, ir <= stepsr/2, ir++,(* Integration over wr *)
    wr = N[(ir-1/2)*2Pi/stepsr,Precision[{a,p,e,\[ScriptCapitalI]}]];
    rp = r[wr];
    Ur = {1,-1,1,-1}*If[wr<Pi,1,-1]*Sqrt[((rp^2+a^2)*En0-a*Lz0)^2-(rp^2-2rp+a^2)*(rp^2+Kc0)];(* Geodesic radial velocity at each quadrant *)
    expr = Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"][wr])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"][wr])+2Pi*n*(ir-1/2)/stepsr)];(* Exponential term with \[CapitalDelta]tr and \[CapitalDelta]\[Phi]r *)
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
      If[ir==1,(* functions of only \[Theta] saved to a list *)
        zp = z[w\[Theta]];
        Uz = {1,1,-1,-1}*If[w\[Theta]<Pi,-1,1]*Sqrt[-((1-zp^2)*a*En0-Lz0)^2+(1-zp^2)*(Kc0-a^2*zp^2)];(* Polar geodesic velocity *)
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
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
        -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
      St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
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
      Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Cdnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + r1p*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - z1p/sin\[Theta]p*vn)*vn;
      Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + r1p*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - z1p/sin\[Theta]p*vn*vmb);
      Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Cdmbmb;
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
        \[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RInrp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRInrp + (Cmmbmb*fmbmb2)*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})];  (* Total over all quadrants *)
      sumMinus += Total[s*\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUprp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUprp + (Cmmbmb*fmbmb2)*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})];
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


(* ::Subsubsection:: *)
(*Corrections*)


TeukolskySpinModeCorrectionNum[l_,m_,n_,k_,orbitCorrection_,\[Delta]\[Omega]_]:=Module[{h1,h2,h3,a,p,e,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,correction,r,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,
    SWSH,SWSHplus,SWSHminus,R,Rplus,Rminus,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    correctionp,lists,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,
    RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,\[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1},
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
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
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
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
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
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RInrp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRInrp + Cmmbmb*fmbmb2*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUprp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUprp + Cmmbmb*fmbmb2*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RInrp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRInrp1 + Cmmbmb*fmbmb2*ddRInrp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RInrp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRInrp + Cmmbmb*fmbmb21*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RUprp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRUprp1 + Cmmbmb*fmbmb2*ddRUprp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RUprp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRUprp + Cmmbmb*fmbmb21*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
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


TeukolskySpinModeCorrection[l_,m_,n_,k_,orbitCorrection_,{angparNew_,RCorrection_}]:=Module[{h1,h2,h3,a,p,e,\[ScriptCapitalI],En,Lz,Kc,En1,Lz1,\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,correction,
    r,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,SWSH,dSWSHd\[Omega],R,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,stepsr,steps\[Theta],\[Theta]list,
    correctionp,lists,ir,i\[Theta],wr,w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,expr,expr1,exp\[Theta],exp\[Theta]1,\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,
    RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,\[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],
    fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1,\[Psi],d\[Lambda]d\[Omega]},
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
  {\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = SetPrecision[orbitCorrection["BLFrequenciesGeo"],30]; (* Geodesic BL frequencies *)
  {\[CapitalOmega]r1,\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Linear corrections to the frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  r = orbitCorrection["TrajectoryGeo"][[2]]; (* Geodesic coordinates r and z=cos(\[Theta]) *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma]  = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]];
  \[Omega]  = m*\[CapitalOmega]\[Phi] + n*\[CapitalOmega]r + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + n*\[CapitalOmega]r1 + k*\[CapitalOmega]\[Theta]1;
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,SetPrecision[a,30],\[Omega],1]; (* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  \[Lambda]1 = d\[Lambda]d\[Omega]*\[Omega]1;
  R = RCorrection[-2,l,m,SetPrecision[a,30],\[Omega],1,\[Lambda],d\[Lambda]d\[Omega],e,p];
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
  stepsr = Max[16*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[Pi  ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[Pi  ]+n)]],
               16*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]tr"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"]'[0   ]+n)]],32];
  steps\[Theta] = Max[ 8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
                8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
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
        d2Sd\[Theta]2 = SWSH[ArcCos[zp]][[3]]; (* second derivative of S from Teukolsky equation *)
        d2S1d\[Theta]2 = dSWSHd\[Omega][ArcCos[zp]][[3]]*\[Omega]1;(* Linear part of second derivative of S from Teukolsky equation *)
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
      Cmnn   = vn^2;
      Cmnmb  = vn*vmb;
      Cmmbmb = vmb^2;
      rho = 1/\[Zeta]; (* Spin coefficients *)
      beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
      pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
      tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
      mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
      gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
      alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
      Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
      Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
      Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
      Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
      Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
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
      Ann0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn + 2*vn1)*vn + Cdnn;
      Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
      AnnrS  = (Srn + rp1*vn)*vn;
      Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
      Anmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
      Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
      AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
      Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
      Ambmb0S  = ((\[CapitalSigma]1/\[CapitalSigma])*vmb + 2*vmb1)*vmb + Cdmbmb;
      Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
      AmbmbrS  = (Srmb + rp1*vmb)*vmb;
      Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RInrp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRInrp + Cmmbmb*fmbmb2*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUprp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUprp + Cmmbmb*fmbmb2*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RInrp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRInrp1 + Cmmbmb*fmbmb2*ddRInrp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RInrp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRInrp + Cmmbmb*fmbmb21*ddRInrp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; (* Total over all quadrants *)
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RUprp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRUprp1 + Cmmbmb*fmbmb2*ddRUprp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RUprp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRUprp + Cmmbmb*fmbmb21*ddRUprp)*expr^{1,-1,1,-1}*exp\[Theta]^{1,1,-1,-1})]; 
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


(* ::Subsection:: *)
(*Spherical*)


TeukolskySpinModeSpherical[l_,m_,k_,orbitCorrection_,s_]:=Module[{h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],correction,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],SWSH,R,\[Lambda],\[ScriptCapitalC]2,rplus,P,\[Epsilon],\[Alpha],W,
    sumPlus,sumMinus,steps\[Theta],correctionp,i\[Theta],w\[Theta],rp,zp,sin\[Theta]p,Ur,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dK,V,dV,RInrp,dRInrp,ddRInrp,dddRInrp,RUprp,dRUprp,ddRUprp,dddRUprp,
    \[Theta]2,S,L2S,L1L2S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus,CMinus},
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
  steps\[Theta] = Max[8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];
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
    Cmnn   = vn^2;
    Cmnmb  = vn*vmb;
    Cmmbmb = vmb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
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
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Cdnn;
    Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Cdmbmb;
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
      \[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RInrp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRInrp + (Cmmbmb*fmbmb2)*ddRInrp)*exp\[Theta]^{1,-1})]; 
    sumMinus += Total[s*\[CapitalSigma]*((
      ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
      ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
      ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
      ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
      ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
      ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
      \[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUprp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUprp + (Cmmbmb*fmbmb2)*ddRUprp)*exp\[Theta]^{1,-1})];
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


TeukolskySpinModeSphericalCorrectionNum[l_,m_,k_,orbitCorrection_,\[Delta]\[Omega]_]:=Module[{h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,correction,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,
    SWSH,SWSHplus,SWSHminus,R,Rplus,Rminus,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,steps\[Theta],correctionp,i\[Theta],w\[Theta],rp,zp,sin\[Theta]p,
    Ur,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,
    \[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,
    fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,
    dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,
    Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,\[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,
    AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1},
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
  steps\[Theta] = Max[8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/4]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/4]+k)]],
               8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];
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
    Cmnn   = vn^2;
    Cmnmb  = vn*vmb;
    Cmmbmb = vmb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
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
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Cdnn;
    Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Cdmbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
    AmbmbrS  = (Srmb + rp1*vmb)*vmb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RInrp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRInrp + Cmmbmb*fmbmb2*ddRInrp)*exp\[Theta]^{1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUprp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUprp + Cmmbmb*fmbmb2*ddRUprp)*exp\[Theta]^{1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RInrp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRInrp1 + Cmmbmb*fmbmb2*ddRInrp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RInrp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRInrp + Cmmbmb*fmbmb21*ddRInrp)*exp\[Theta]^{1,-1})]; 
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RUprp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRUprp1 + Cmmbmb*fmbmb2*ddRUprp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RUprp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRUprp + Cmmbmb*fmbmb21*ddRUprp)*exp\[Theta]^{1,-1})]; 
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
    "steps\[Theta]"->steps\[Theta]
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]


(*Options[TeukolskySpinModeSphericalCorrection] = {WorkingPrecision->30};

TeukolskySpinModeSphericalCorrection[l_,m_,k_,orbitCorrection_,{angparNew_,TeukolskySolverHS1spin_},OptionsPattern[]]:=Module[{h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,
    correction,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,SWSH,dSWSHd\[Omega],R,\[Lambda],d\[Lambda]d\[Omega],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,steps\[Theta],correctionp,i\[Theta],w\[Theta],rp,zp,sin\[Theta]p,
    Ur,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,
    \[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,
    fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,
    Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,
    \[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1,\[Omega]prec,\[Theta]p},
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_] := 1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+
               a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  Print["Calculating l = "<>ToString[l]<>", m = "<>ToString[m]<>", k = "<>ToString[k]<>" mode"];
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  x = orbitCorrection["Inclination"];
  En = orbitCorrection["Ehat"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzhat"];
  Kc = orbitCorrection["Khat"];
  En1 = orbitCorrection["ES"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["LS"];
  {\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Coordinate frequencies *)
  {\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Coordinate frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]]; 
  \[Omega] = m*\[CapitalOmega]\[Phi] + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  \[Omega]prec = m*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"]
        + k*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)"];
  (*\[Omega]prec = SetPrecision[\[Omega],40];*)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + k*\[CapitalOmega]\[Theta]1;
  (*Print["Calculating angular functions"<>ToString@AbsoluteTiming[{\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1];][[1]]];*)(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  {\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega]}=angparNew[-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1];(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  \[Lambda]1 = d\[Lambda]d\[Omega]*\[Omega]1;
  Print["Calculating radial functions"<>ToString@AbsoluteTiming[R = TeukolskySolverHS1spin[p,-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1,\[Lambda],d\[Lambda]d\[Omega]]][[1]]];
  (*R = TeukolskySolverHS1spin[p,-2,l,m,SetPrecision[a,OptionValue[WorkingPrecision]+5],\[Omega]prec,1,\[Lambda],d\[Lambda]d\[Omega]];*)
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
  steps\[Theta] = Max[8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/2]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/2]+k)]],
               8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];
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
  (*{{{RInrp,dRInrp},{RInrp1,dRInrp1}},{{RUprp,dRUprp},{RUprp1,dRUprp1}}} = R*{{{1,1},{\[Omega]1,\[Omega]1}},{{1,1},{\[Omega]1,\[Omega]1}}};*)
  {RInrp,dRInrp} = R[[1,1]];
  (*RInrp    = R[[1,1,1]]; (* radial function *)
  dRInrp   = R[[1,1,2]];*)
  ddRInrp  = (V*RInrp + d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  {RInrp1,dRInrp1} = R[[1,2]]*\[Omega]1;
  (*RInrp1   = R[[1,2,1]]*\[Omega]1; (* Linear part of the radial function *)
  dRInrp1  = R[[1,2,2]]*\[Omega]1;*)
  ddRInrp1 = (V*RInrp1 + V1*RInrp + d\[CapitalDelta]*dRInrp1)/\[CapitalDelta];
  dddRInrp = 1/\[CapitalDelta] (dV*RInrp + (V + 2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
  {RUprp,dRUprp} = R[[2,1]];
  (*RUprp    = R[[2,1,1]];
  dRUprp   = R[[2,1,2]];*)
  ddRUprp  = (V*RUprp + d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
  {RUprp1,dRUprp1} = R[[2,2]]*\[Omega]1;
  (*RUprp1   = R[[2,2,1]]*\[Omega]1;
  dRUprp1  = R[[2,2,2]]*\[Omega]1;*)
  ddRUprp1 = (V*RUprp1 + V1*RUprp + d\[CapitalDelta]*dRUprp1)/\[CapitalDelta];
  dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
  Do[(* integration over w\[Theta] *)
    w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]];
    zp=z[w\[Theta]];(* functions of only w\[Theta] saved to lists *)
    Uz={1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
    exp\[Theta]=Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
    sin\[Theta]p=Sqrt[1-zp^2];
    \[Theta]p=ArcCos[zp];
    {S,dSd\[Theta],d2Sd\[Theta]2}=SWSH[\[Theta]p];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
    {S1,dS1d\[Theta],d2S1d\[Theta]2}=dSWSHd\[Omega][\[Theta]p]*\[Omega]1;  (*  Linear part of S(\[Theta](z))  *)
    d3Sd\[Theta]3=-(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
           -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
    L2S=dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
    L2S1=dS1d\[Theta]-(S1 (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2) \[Omega]1))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
    dL2Sd\[Theta]=d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
    L1L2S=d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
    L1L2S1=d2S1d\[Theta]2+(dS1d\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2) \[Omega]1))/sin\[Theta]p + 
           S1 (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) \[Omega]1 - a^2 (-1+zp^2) 2*\[Omega]*\[Omega]1);  
    dL1L2Sd\[Theta]=d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
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
    Cmnn   = vn^2;
    Cmnmb  = vn*vmb;
    Cmmbmb = vmb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
      -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
    St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
    S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
    rp1 = {correctionp[[i\[Theta]]]["rS"],correctionp[[-i\[Theta]]]["rS"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
    zp1 = {correctionp[[i\[Theta]]]["zS"],-correctionp[[-i\[Theta]]]["zS"]};
    Urp1 = {correctionp[[i\[Theta]]]["UrS"],correctionp[[-i\[Theta]]]["UrS"]};
    Uzp1 = {correctionp[[i\[Theta]]]["UzS"],-correctionp[[-i\[Theta]]]["UzS"]};
    \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
    exp1 = I*{(\[Omega]*correctionp[[ i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ i\[Theta]]]["\[CapitalDelta]\[Phi]S"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
              (\[Omega]*correctionp[[-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};
    vn1  = -( ((2*rp*rp1)*En) - ((rp^2+a^2)*En - a*Lz + Ur)/(\[CapitalSigma])*\[CapitalSigma]1 + 
             ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]) + Urp1))/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
    vmb1 = ( (-I*(-2*a*zp*zp1*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp])) + Uzp1))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Cdnn;
    Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Cdmbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
    AmbmbrS  = (Srmb + rp1*vmb)*vmb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RInrp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRInrp + Cmmbmb*fmbmb2*ddRInrp)*exp\[Theta]^{1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUprp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUprp + Cmmbmb*fmbmb2*ddRUprp)*exp\[Theta]^{1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RInrp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRInrp1 + Cmmbmb*fmbmb2*ddRInrp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RInrp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRInrp + Cmmbmb*fmbmb21*ddRInrp)*exp\[Theta]^{1,-1})]; 
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RUprp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRUprp1 + Cmmbmb*fmbmb2*ddRUprp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RUprp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRUprp + Cmmbmb*fmbmb21*ddRUprp)*exp\[Theta]^{1,-1})];,
    {i\[Theta], 1, steps\[Theta]/2}
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
    "S"->N[SWSH[Pi/2,0][[1]]],
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
    "steps\[Theta]"->steps\[Theta]
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]*)


Options[TeukolskySpinModeSphericalCorrection] = {WorkingPrecision->30};

TeukolskySpinModeSphericalCorrection[l_,m_,k_,orbitCorrection_,orbitDerivatives_,{angparNew_,TeukolskySolverHS1spin_},OptionsPattern[]]:=Module[{
    h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,dEndr,dLzdr,dEndx,dLzdx,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,d\[CapitalOmega]\[Theta]dr,d\[CapitalOmega]\[Phi]dr,d\[CapitalOmega]\[Theta]dx,d\[CapitalOmega]\[Phi]dx,correction,derivatives,z,\[CapitalGamma],\[CapitalGamma]1,d\[CapitalGamma]dr,d\[CapitalGamma]dx,
    \[Omega],\[Omega]prec,\[Omega]1,d\[Omega]dr,d\[Omega]dx,\[Lambda],d\[Lambda]d\[Omega],SWSH,dSWSHd\[Omega],R,\[ScriptCapitalC]2,d\[ScriptCapitalC]2d\[Omega],rplus,P,\[Epsilon],\[Alpha],d\[Alpha]d\[Omega],
    steps\[Theta],correctionp,derivativesp,sumPlus0,sumMinus0,sumPlus1,sumMinus1,dsumPlusdr,dsumMinusdr,dsumPlusdx,dsumMinusdx,i\[Theta],w\[Theta],
    rp,zp,\[Theta]p,sin\[Theta]p,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,dKd\[Omega],dKdr,d2Kdrd\[Omega],V,dVd\[Omega],dVdr,
    RIn,dRIndr,d2RIndr2,d3RIndr3,dRInd\[Omega],d2RIndrd\[Omega],d3RIndr2d\[Omega],RUp,dRUpdr,d2RUpdr2,d3RUpdr3,dRUpd\[Omega],d2RUpdrd\[Omega],d3RUpdr2d\[Omega],
    \[Theta]2,S,dSd\[Theta],d2Sd\[Theta]2,d3Sd\[Theta]3,dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega],L2S,dL2Sd\[Theta],dL2Sd\[Omega],L1L2S,dL1L2Sd\[Theta],dL1L2Sd\[Omega],
    \[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,dfnn0d\[Omega],dfnmb0d\[Omega],dfnmb1d\[Omega],dfmbmb0d\[Omega],dfmbmb1d\[Omega],dfmbmb2d\[Omega],dfnn0dr,
    dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],
    ul,un,um,umb,Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,
    St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,
    rp1,zp1,dzpdr,d\[Theta]pdr,dzpdx,d\[Theta]pdx,Urp1,Uzp1,dUzpdr,dUzpdx,\[CapitalSigma]1,d\[CapitalSigma]dr,d\[CapitalSigma]dx,exp1,dexpdr,dexpdx,un1,dundr,dundx,umb1,dumbdr,dumbdx,
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
  \[Omega]prec = m*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)"]
        + k*KerrGeodesics`OrbitalFrequencies`KerrGeoFrequencies[SetPrecision[a,OptionValue[WorkingPrecision]+5],SetPrecision[p,OptionValue[WorkingPrecision]+5],0,SetPrecision[x,OptionValue[WorkingPrecision]+5]]["\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)"]; (* Frequency with higher precision *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + k*\[CapitalOmega]\[Theta]1;
  d\[Omega]dr = m*d\[CapitalOmega]\[Phi]dr + k*d\[CapitalOmega]\[Theta]dr;
  d\[Omega]dx = m*d\[CapitalOmega]\[Phi]dx + k*d\[CapitalOmega]\[Theta]dx;
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
  (* number of steps for w\[Theta] integration *)
  steps\[Theta] = Max[8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/2]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/2]+k)]],
               8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];
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
    {S,dSd\[Theta],d2Sd\[Theta]2}=SWSH[\[Theta]p];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
    {dSd\[Omega],d2Sd\[Theta]d\[Omega],d3Sd\[Theta]2d\[Omega]}=dSWSHd\[Omega][\[Theta]p];  (*  \[Omega]-derivative of S(\[Theta](z))  *)
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
    Cmnn   = un^2;
    Cmnmb  = un*umb;
    Cmmbmb = umb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*un-Sln*2*Re[gamma]*un-2*Re[Snmb*((Conjugate[alpha]+beta)*un-mu*um)]);
    Cdmbmb= (Scd\[Gamma]mbdc*umb-Snmb*(-pi*ul)-Slmb*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)+Smmb*(-(-alpha+Conjugate[beta])*umb));
    Cdnmb = (Scd\[Gamma]ndc*umb+Scd\[Gamma]mbdc*un-Sln*(Conjugate[tau]*un-(Conjugate[gamma]-gamma)*umb)-Snmb*(Conjugate[rho]*un-mu*ul-(Conjugate[alpha]-beta)*umb)
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
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*un + 2*un1)*un + Cdnn; (* Source term *)
    Annt\[Phi]S = (St\[Phi]n + exp1*un)*un;
    AnnrS  = (Srn + rp1*un)*un;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*un)*un;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*un*umb + un1*umb + un*umb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*umb + St\[Phi]mb*un)/2 + exp1*un*umb);
    AnmbrS  = ((Srn*umb + Srmb*un)/2 + rp1*un*umb);
    Anmb\[Theta]S  = ((S\[Theta]n*umb + S\[Theta]mb*un)/2 - zp1/sin\[Theta]p*un*umb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*umb + 2*umb1)*umb + Cdmbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*umb)*umb;
    AmbmbrS  = (Srmb + rp1*umb)*umb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*umb)*umb;
      {Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RIn - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRIndr + Cmmbmb*fmbmb2*d2RIndr2)*exp\[Theta]^{1,-1})], (* Totral of all quadrants *)
       Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUpdr + Cmmbmb*fmbmb2*d2RUpdr2)*exp\[Theta]^{1,-1})], 
       Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RIn*fnn0 + AnnrS*(dRIndr*fnn0 + RIn*dfnn0dr) + Ann\[Theta]S*RIn*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RIn*fnmb0 + AnmbrS*( RIn*dfnmb0dr +  dRIndr*fnmb0) + Anmb\[Theta]S* RIn*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRIndr*fnmb1 + AnmbrS*(dRIndr*dfnmb1dr + d2RIndr2*fnmb1) + Anmb\[Theta]S*dRIndr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RIn*fmbmb0 + AmbmbrS*(  dRIndr*fmbmb0 +   RIn*dfmbmb0dr) + Ambmb\[Theta]S*  RIn*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRIndr*fmbmb1 + AmbmbrS*( d2RIndr2*fmbmb1 +  dRIndr*dfmbmb1dr) + Ambmb\[Theta]S* dRIndr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RIndr2*fmbmb2 + AmbmbrS*(d3RIndr3*fmbmb2 + d2RIndr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RIndr2*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*dRInd\[Omega] - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*d2RIndrd\[Omega] + Cmmbmb*fmbmb2*d3RIndr2d\[Omega]
          + (Cmnn*dfnn0d\[Omega] + Cmnmb*dfnmb0d\[Omega] + Cmmbmb*dfmbmb0d\[Omega])*RIn - (Cmnmb*dfnmb1d\[Omega] + Cmmbmb*dfmbmb1d\[Omega])*dRIndr + Cmmbmb*dfmbmb2d\[Omega]*d2RIndr2)*\[Omega]1*exp\[Theta]^{1,-1})], 
       Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUp*fnn0 + AnnrS*(dRUpdr*fnn0 + RUp*dfnn0dr) + Ann\[Theta]S*RUp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUp*fnmb0 + AnmbrS*( RUp*dfnmb0dr +  dRUpdr*fnmb0) + Anmb\[Theta]S* RUp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUpdr*fnmb1 + AnmbrS*(dRUpdr*dfnmb1dr + d2RUpdr2*fnmb1) + Anmb\[Theta]S*dRUpdr*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUp*fmbmb0 + AmbmbrS*(  dRUpdr*fmbmb0 +   RUp*dfmbmb0dr) + Ambmb\[Theta]S*  RUp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUpdr*fmbmb1 + AmbmbrS*( d2RUpdr2*fmbmb1 +  dRUpdr*dfmbmb1dr) + Ambmb\[Theta]S* dRUpdr*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*d2RUpdr2*fmbmb2 + AmbmbrS*(d3RUpdr3*fmbmb2 + d2RUpdr2*dfmbmb2dr) + Ambmb\[Theta]S*d2RUpdr2*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*dRUpd\[Omega] - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*d2RUpdrd\[Omega] + Cmmbmb*fmbmb2*d3RUpdr2d\[Omega]
          + (Cmnn*dfnn0d\[Omega] + Cmnmb*dfnmb0d\[Omega] + Cmmbmb*dfmbmb0d\[Omega])*RUp - (Cmnmb*dfnmb1d\[Omega] + Cmmbmb*dfmbmb1d\[Omega])*dRUpdr + Cmmbmb*dfmbmb2d\[Omega]*d2RUpdr2)*\[Omega]1*exp\[Theta]^{1,-1})], 
      Total[\[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Cmnn*fnn0*RIn + Cmnmb*(fnmb0*RIn - fnmb1*dRIndr) + Cmmbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundr*fnn0*RIn + (un*dumbdr+dundr*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdr*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Cmnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RIn + fnn0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)) + 
               Cmnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RIn + fnmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRIndr - fnmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr)) + 
              Cmmbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RIn + fmbmb0*(dRIndr + dRInd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRIndr - fmbmb1*(d2RIndr2 + d2RIndrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RIndr2 + fmbmb2*(d3RIndr3 + d3RIndr2d\[Omega]*d\[Omega]dr)))) 
              )*exp\[Theta]^{1,-1}], 
       Total[\[CapitalSigma]*((d\[CapitalSigma]dr/\[CapitalSigma] + dexpdr)*((Cmnn*fnn0*RUp + Cmnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Cmmbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundr*fnn0*RUp + (un*dumbdr+dundr*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdr*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Cmnn*((dfnn0dr + dfnn0d\[Theta]*d\[Theta]pdr + dfnn0d\[Omega]*d\[Omega]dr)*RUp + fnn0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)) + 
               Cmnmb*((dfnmb0dr + dfnmb0d\[Theta]*d\[Theta]pdr + dfnmb0d\[Omega]*d\[Omega]dr)*RUp + fnmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfnmb1dr + dfnmb1d\[Theta]*d\[Theta]pdr + dfnmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fnmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr)) + 
              Cmmbmb*((dfmbmb0dr + dfmbmb0d\[Theta]*d\[Theta]pdr + dfmbmb0d\[Omega]*d\[Omega]dr)*RUp + fmbmb0*(dRUpdr + dRUpd\[Omega]*d\[Omega]dr)
                    - (dfmbmb1dr + dfmbmb1d\[Theta]*d\[Theta]pdr + dfmbmb1d\[Omega]*d\[Omega]dr)*dRUpdr - fmbmb1*(d2RUpdr2 + d2RUpdrd\[Omega]*d\[Omega]dr) + 
                      (dfmbmb2dr + dfmbmb2d\[Theta]*d\[Theta]pdr + dfmbmb2d\[Omega]*d\[Omega]dr)*d2RUpdr2 + fmbmb2*(d3RUpdr3 + d3RUpdr2d\[Omega]*d\[Omega]dr)))) 
              )*exp\[Theta]^{1,-1}], 
      Total[\[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Cmnn*fnn0*RIn + Cmnmb*(fnmb0*RIn - fnmb1*dRIndr) + Cmmbmb*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) +
              ((2*un*dundx*fnn0*RIn + (un*dumbdx+dundx*umb)*(fnmb0*RIn - fnmb1*dRIndr) + 2*umb*dumbdx*(fmbmb0*RIn - fmbmb1*dRIndr + fmbmb2*d2RIndr2))) + 
              ((Cmnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RIn + fnn0*(dRInd\[Omega]*d\[Omega]dx)) + 
               Cmnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RIn + fnmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRIndr - fnmb1*(d2RIndrd\[Omega]*d\[Omega]dx)) + 
              Cmmbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RIn + fmbmb0*(dRInd\[Omega]*d\[Omega]dx)
                    - (dfmbmb1d\[Theta]*d\[Theta]pdx + dfmbmb1d\[Omega]*d\[Omega]dx)*dRIndr - fmbmb1*(d2RIndrd\[Omega]*d\[Omega]dx) + 
                      (dfmbmb2d\[Theta]*d\[Theta]pdx + dfmbmb2d\[Omega]*d\[Omega]dx)*d2RIndr2 + fmbmb2*(d3RIndr2d\[Omega]*d\[Omega]dx)))) 
              )*exp\[Theta]^{1,-1}], 
       Total[\[CapitalSigma]*((d\[CapitalSigma]dx/\[CapitalSigma] + dexpdx)*((Cmnn*fnn0*RUp + Cmnmb*(fnmb0*RUp - fnmb1*dRUpdr) + Cmmbmb*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) +
              ((2*un*dundx*fnn0*RUp + (un*dumbdx+dundx*umb)*(fnmb0*RUp - fnmb1*dRUpdr) + 2*umb*dumbdx*(fmbmb0*RUp - fmbmb1*dRUpdr + fmbmb2*d2RUpdr2))) + 
              ((Cmnn*((dfnn0d\[Theta]*d\[Theta]pdx + dfnn0d\[Omega]*d\[Omega]dx)*RUp + fnn0*(dRUpd\[Omega]*d\[Omega]dx)) + 
               Cmnmb*((dfnmb0d\[Theta]*d\[Theta]pdx + dfnmb0d\[Omega]*d\[Omega]dx)*RUp + fnmb0*(dRUpd\[Omega]*d\[Omega]dx)
                    - (dfnmb1d\[Theta]*d\[Theta]pdx + dfnmb1d\[Omega]*d\[Omega]dx)*dRUpdr - fnmb1*(d2RUpdrd\[Omega]*d\[Omega]dx)) + 
              Cmmbmb*((dfmbmb0d\[Theta]*d\[Theta]pdx + dfmbmb0d\[Omega]*d\[Omega]dx)*RUp + fmbmb0*(dRUpd\[Omega]*d\[Omega]dx)
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
  <|
    "l" -> l,
    "m" -> m,
    "k" -> k,
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
    |>,
    "steps\[Theta]" -> steps\[Theta]
  |>
]


(*TeukolskySpinModeSphericalCorrection2[l_,m_,k_,orbitCorrection_,{angparNew_,TeukolskySolverHS1spin_}]:=Module[{h1,h2,a,p,e,x,En,Lz,Kc,En1,Lz1,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi],\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1,
    correction,z,\[CapitalGamma],\[CapitalGamma]1,\[Omega],\[Omega]1,SWSH,SWSHS,R,\[Lambda],\[Lambda]1,\[ScriptCapitalC]2,\[ScriptCapitalC]21,rplus,P,\[Epsilon],\[Alpha],\[Alpha]1,W,W1,sumPlus0,sumMinus0,sumPlus1,sumMinus1,steps\[Theta],correctionp,i\[Theta],w\[Theta],rp,zp,sin\[Theta]p,
    Ur,Uz,exp\[Theta],\[CapitalDelta],d\[CapitalDelta],K,K1,dK,dK1,V,V1,dV,RInrp,dRInrp,ddRInrp,RInrp1,dRInrp1,ddRInrp1,dddRInrp,RUprp,dRUprp,ddRUprp,RUprp1,dRUprp1,ddRUprp1,dddRUprp,
    \[Theta]2,S,S1,L2S,L2S1,L1L2S,L1L2S1,dSd\[Theta],dS1d\[Theta],d2Sd\[Theta]2,d2S1d\[Theta]2,d3Sd\[Theta]3,dL2Sd\[Theta],dL1L2Sd\[Theta],\[Zeta],\[Zeta]bar,\[CapitalSigma],fnn0,fnmb0,fnmb1,fmbmb0,fmbmb1,fmbmb2,fnn01,fnmb01,fnmb11,
    fmbmb01,fmbmb11,fmbmb21,dfnn0dr,dfnmb0dr,dfnmb1dr,dfmbmb0dr,dfmbmb1dr,dfmbmb2dr,dfnn0d\[Theta],dfnmb0d\[Theta],dfnmb1d\[Theta],dfmbmb0d\[Theta],dfmbmb1d\[Theta],dfmbmb2d\[Theta],vl,vn,vm,vmb,
    Sln,Slmb,Snm,Snmb,Smmb,Cmnn,Cmnmb,Cmmbmb,rho,beta,pi,alpha,mu,gamma,tau,Scd\[Gamma]ndc,Scd\[Gamma]mbdc,Cdnn,Cdnmb,Cdmbmb,St\[Phi]n,St\[Phi]mb,Srn,Srmb,S\[Theta]n,S\[Theta]mb,rp1,zp1,Urp1,Uzp1,
    \[CapitalSigma]1,exp1,vn1,vmb1,Ann0S,Annt\[Phi]S,AnnrS,Ann\[Theta]S,Anmb0S,Anmbt\[Phi]S,AnmbrS,Anmb\[Theta]S,Ambmb0S,Ambmbt\[Phi]S,AmbmbrS,Ambmb\[Theta]S,CPlus0,CMinus0,CPlus1,CMinus1},
  h1[r_,z_] := (r (-3 a^2 r^2 z^2+a^4 z^4+Kc (r^2-3 a^2 z^2)))/(Sqrt[Kc] (r^2+a^2 z^2)^3);
  h2[r_,z_] := 1/(Sqrt[Kc] (r^2+a^2 z^2)^3) (-En Lz r^6+a^4 En Lz r^2 z^4+a^2 En Lz r^4 (-2+z^2)-a^6 En Lz z^4 (-2+z^2)+
               a^7 En^2 z^4 (-1+z^2)+a r^3 (Lz^2 r+Kc (-1+z^2)+Kc r (-1+2 z^2)+r^3 (z^2-En^2 (-1+z^2)))+a^3 (-En^2 r^4 (-1+z^2)+
               r z^2 (2 r^3+2 Kc r z^2-3 Kc (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Kc-Lz^2+r ( (-1+z^2)+r (2-z^2+En^2 (-1+z^2)))));
  a = orbitCorrection["a"];(* Orbital parameters *)
  p = orbitCorrection["p"];
  e = orbitCorrection["e"];
  x = orbitCorrection["Inclination"];
  En = orbitCorrection["Ehat"]; (* Geodesic constants of motion *)
  Lz = orbitCorrection["Lzhat"];
  Kc = orbitCorrection["Khat"];
  En1 = orbitCorrection["ES"]; (* Linear corrections to the constants of motion *)
  Lz1 = orbitCorrection["LS"];
  {\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]} = orbitCorrection["BLFrequenciesGeo"]; (* Coordinate frequencies *)
  {\[CapitalOmega]\[Theta]1,\[CapitalOmega]\[Phi]1} = orbitCorrection["BLFrequenciesCorrection"]; (* Coordinate frequencies *)
  correction = orbitCorrection["OrbitCorrection"]; (* function containing corrections to the trajectory *)
  z[wz_] := Cos[orbitCorrection["TrajectoryGeo"][[3]][wz]];
  \[CapitalGamma] = orbitCorrection["MinoFrequenciesGeo"][[1]]; (* Geodesic average rate of change of BL time in Mino time and the linear correction *)
  \[CapitalGamma]1 = orbitCorrection["MinoFrequenciesCorrection"][[1]]; 
  \[Omega] = m*\[CapitalOmega]\[Phi] + k*\[CapitalOmega]\[Theta]; (* Geodesic frequency and the linear correction *)
  \[Omega]1 = m*\[CapitalOmega]\[Phi]1 + k*\[CapitalOmega]\[Theta]1;
  {\[Lambda],\[Lambda]1,SWSH,SWSHS}=angparNew[-2,l,m,SetPrecision[a,48],SetPrecision[\[Omega],48],SetPrecision[\[Omega]1,48]];(* Polar and radial functions and the eigenvalue for geodesic frequency and linear corrections *)
  R = TeukolskySolverHS1spin[p,-2,l,m,SetPrecision[a,48],SetPrecision[\[Omega],48],SetPrecision[\[Omega]1,48],\[Lambda],\[Lambda]1];
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
  steps\[Theta] = Max[8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[Pi/2]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[Pi/2]+k)]],
               8*Ceiling[Abs[(\[Omega]*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"]'[0   ]-m*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]'[0   ]+k)]],32];
  Print[ToString[steps\[Theta]]<>" steps in w\[Theta]"];
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
  RInrp    = R[[1,1,1]]; (* radial function *)
  dRInrp   = R[[1,1,2]];
  ddRInrp  = (V*RInrp + d\[CapitalDelta]*dRInrp)/\[CapitalDelta];  (* second derivative of radial function from Teukolsky equation *)
  RInrp1   = R[[1,2,1]]; (* Linear part of the radial function *)
  dRInrp1  = R[[1,2,2]];
  ddRInrp1 = (V*RInrp1 + V1*RInrp + d\[CapitalDelta]*dRInrp1)/\[CapitalDelta];
  dddRInrp = 1/\[CapitalDelta] (dV*RInrp + (V + 2)*dRInrp);  (* third derivative of radial function from Teukolsky equation *)
  RUprp    = R[[2,1,1]];
  dRUprp   = R[[2,1,2]];
  ddRUprp  = (V*RUprp + d\[CapitalDelta]*dRUprp)/\[CapitalDelta];  
  RUprp1   = R[[2,2,1]];
  dRUprp1  = R[[2,2,2]];
  ddRUprp1 = (V*RUprp1 + V1*RUprp + d\[CapitalDelta]*dRUprp1)/\[CapitalDelta];
  dddRUprp = 1/\[CapitalDelta] (dV*RUprp+(V+2)*dRUprp);
  For[i\[Theta] = 1, i\[Theta] <= steps\[Theta]/2, i\[Theta]++,(* integration over w\[Theta] *)
    w\[Theta]=N[(i\[Theta]-1/2)*2Pi/steps\[Theta],Precision[{a,p,e,x}]];
    zp=z[w\[Theta]];(* functions of only w\[Theta] saved to lists *)
    Uz={1,-1}*(-1)*Sqrt[-((1-zp^2)*a*En-Lz)^2+(1-zp^2)*(Kc-a^2*zp^2)];(* Polar geodesic velocity *)
    exp\[Theta]=Exp[I*(\[Omega]*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]])-m*(orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"][w\[Theta]])+2Pi*k*(i\[Theta]-1/2)/steps\[Theta])];
    sin\[Theta]p=Sqrt[1-zp^2];
    S=SWSH[ArcCos[zp]][[1]];  (*  Spin-weighted spheroidal harmonics S(\[Theta](z))  *)
    S1=SWSHS[ArcCos[zp]][[1]];  (*  Linear part of S(\[Theta](z))  *)
    dSd\[Theta]=SWSH[ArcCos[zp]][[2]]; (* First derivative of S wrt \[Theta] *)
    dS1d\[Theta]=SWSHS[ArcCos[zp]][[2]]; (* Linear part of first derivative of S wrt \[Theta] *)
    d2Sd\[Theta]2=SWSH[ArcCos[zp]][[3]]; (* second derivative of S from Teukolsky equation *)
    d2S1d\[Theta]2=SWSHS[ArcCos[zp]][[3]];(* Linear part of second derivative of S from Teukolsky equation *)
    d3Sd\[Theta]3=-(1/sin\[Theta]p^3)2 (-2+m zp+a zp (-1+zp^2) \[Omega]) (m+a \[Omega]-zp (2+a zp \[Omega]))*S
           -(-(a*\[Omega])^2*(1-zp^2)-(m-2*zp)^2/(1-zp^2)+4a*\[Omega]*zp-2+2*m*a*\[Omega]+\[Lambda]-1/(1-zp^2))*dSd\[Theta]-zp/sin\[Theta]p*d2Sd\[Theta]2; (*third derivative from derivative of second derivative*)
    L2S=dSd\[Theta]-(S (m-2 zp+a (-1+zp^2) \[Omega]))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
    L2S1=dS1d\[Theta]-(S1 (m-2 zp+a (-1+zp^2) \[Omega])+S (a (-1+zp^2) \[Omega]1))/sin\[Theta]p;(* Operators acting on S(\[Theta]) *)
    dL2Sd\[Theta]=d2Sd\[Theta]2+1/(-1+zp^2) (dSd\[Theta] sin\[Theta]p (m-2 zp+a (-1+zp^2) \[Omega])+S (2-m zp+a zp (-1+zp^2) \[Omega]));
    L1L2S=d2Sd\[Theta]2+(dSd\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega]))/sin\[Theta]p+S (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2);  
    L1L2S1=d2S1d\[Theta]2+(dS1d\[Theta] (-2 m+3 zp-2 a (-1+zp^2) \[Omega])+dSd\[Theta] (-2 a (-1+zp^2) \[Omega]1))/sin\[Theta]p + 
           S1 (-2-(m (m-2 zp))/(-1+zp^2)-2 a (m-2 zp) \[Omega]-a^2 (-1+zp^2) \[Omega]^2) + S (-2 a (m-2 zp) \[Omega]1 - a^2 (-1+zp^2) 2*\[Omega]*\[Omega]1);  
    dL1L2Sd\[Theta]=d3Sd\[Theta]3+1/(-1+zp^2) dSd\[Theta] (5-m^2-2 zp^2-2 a (m-3 zp) (-1+zp^2) \[Omega]-a^2 (-1+zp^2)^2 \[Omega]^2)+
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
    Cmnn   = vn^2;
    Cmnmb  = vn*vmb;
    Cmmbmb = vmb^2;
    rho = 1/\[Zeta]; (* Spin coefficients *)
    beta = -(zp/(2*\[Zeta]bar Sqrt[2]*sin\[Theta]p));
    pi = -((I a sin\[Theta]p)/(\[Zeta]^2 Sqrt[2]));
    tau = (I a sin\[Theta]p)/(Sqrt[2] \[CapitalSigma]);
    mu = \[CapitalDelta]/(2 \[Zeta]^2 \[Zeta]bar);
    gamma = (a^2-rp+I a (-1+rp) zp)/(2 \[Zeta]^2 \[Zeta]bar);
    alpha = -((-rp zp-I a (-2+zp^2))/(2 \[Zeta]^2 Sqrt[2]sin\[Theta]p));
    Scd\[Gamma]ndc = -Sln*2*Re[gamma](*-2*Re[Snmb*(-Conjugate[pi]+Conjugate[alpha]+beta)]*)-Smmb*(-mu+Conjugate[mu]);
    Scd\[Gamma]mbdc = -Sln*(pi+Conjugate[tau])-Snmb*Conjugate[rho]-Slmb*(-Conjugate[gamma]+gamma-mu)-Smmb*(-alpha+Conjugate[beta]);
    Cdnn  = (Scd\[Gamma]ndc*vn-Sln*2*Re[gamma]*vn-2*Re[Snmb*((Conjugate[alpha]+beta)*vn-mu*vm)]);
    Cdmbmb= (Scd\[Gamma]mbdc*vmb-Snmb*(-pi*vl)-Slmb*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)+Smmb*(-(-alpha+Conjugate[beta])*vmb));
    Cdnmb = (Scd\[Gamma]ndc*vmb+Scd\[Gamma]mbdc*vn-Sln*(Conjugate[tau]*vn-(Conjugate[gamma]-gamma)*vmb)-Snmb*(Conjugate[rho]*vn-mu*vl-(Conjugate[alpha]-beta)*vmb)
      -Snm*(-(-alpha+Conjugate[beta])*vmb)-Snmb*(-Conjugate[pi]*vmb-pi*vm)-Slmb*(2*Re[gamma]*vn)+Smmb*((alpha+Conjugate[beta])*vn-Conjugate[mu]*vmb))/2;
    St\[Phi]n  = -I*K/(2\[CapitalSigma])*Sln+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[CapitalSigma])*(\[Zeta]*Snmb-\[Zeta]bar*Snm);
    St\[Phi]mb = -I*K*(1/\[CapitalDelta]*Snmb+1/(2\[CapitalSigma])*Slmb)+(a*\[Omega]*sin\[Theta]p-m/sin\[Theta]p)/(Sqrt[2]*\[Zeta])*Smmb;
    Srn  = \[CapitalDelta]/(2\[CapitalSigma])*Sln;
    Srmb = -Snmb+\[CapitalDelta]/(2\[CapitalSigma])*Slmb;
    S\[Theta]n  = -(Snmb*\[Zeta]+Snm*\[Zeta]bar)/(Sqrt[2]*\[CapitalSigma]);
    S\[Theta]mb = Smmb/(Sqrt[2]*\[Zeta]);
    rp1 = {correctionp[[i\[Theta]]]["rS"],correctionp[[-i\[Theta]]]["rS"]}; (* Corrections to the coordinates and four-velocity for each quadrant *)
    zp1 = {correctionp[[i\[Theta]]]["zS"],-correctionp[[-i\[Theta]]]["zS"]};
    Urp1 = {correctionp[[i\[Theta]]]["UrS"],correctionp[[-i\[Theta]]]["UrS"]};
    Uzp1 = {correctionp[[i\[Theta]]]["UzS"],-correctionp[[-i\[Theta]]]["UzS"]};
    \[CapitalSigma]1 = 2*(rp*rp1+a^2*zp*zp1); (* Linear correction to \[CapitalSigma] *)
    exp1 = I*{(\[Omega]*correctionp[[ i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[ i\[Theta]]]["\[CapitalDelta]\[Phi]S"])+\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]],
              (\[Omega]*correctionp[[-i\[Theta]]]["\[CapitalDelta]tS"]-m*correctionp[[-i\[Theta]]]["\[CapitalDelta]\[Phi]S"])-\[Omega]1*orbitCorrection["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][w\[Theta]]};
    vn1  = -( ((2*rp*rp1)*En) - ((rp^2+a^2)*En - a*Lz + Ur)/(\[CapitalSigma])*\[CapitalSigma]1 + 
             ((rp^2+a^2)*(En1-h1[rp,zp]) - a*(Lz1+h2[rp,zp]) + Urp1))/(2*\[CapitalSigma]);(* Linear parts of the four-velocity in Kinnersley tetrad *)
    vmb1 = ( (-I*(-2*a*zp*zp1*En)) - (-I*(a*sin\[Theta]p^2*En - Lz) + Uz)*(-zp*zp1/sin\[Theta]p^2 + (rp1-I*a*zp1)/\[Zeta]) + 
             (-I*(a*sin\[Theta]p^2*(En1-h1[rp,zp]) - (Lz1+h2[rp,zp])) + Uzp1))/(-Sqrt[2]*sin\[Theta]p*\[Zeta]);
    Ann0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn + 2*vn1)*vn + Cdnn;
    Annt\[Phi]S = (St\[Phi]n + exp1*vn)*vn;
    AnnrS  = (Srn + rp1*vn)*vn;
    Ann\[Theta]S  = (S\[Theta]n - zp1/sin\[Theta]p*vn)*vn;
    Anmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vn*vmb + vn1*vmb + vn*vmb1 + Cdnmb);
    Anmbt\[Phi]S = ((St\[Phi]n*vmb + St\[Phi]mb*vn)/2 + exp1*vn*vmb);
    AnmbrS  = ((Srn*vmb + Srmb*vn)/2 + rp1*vn*vmb);
    Anmb\[Theta]S  = ((S\[Theta]n*vmb + S\[Theta]mb*vn)/2 - zp1/sin\[Theta]p*vn*vmb);
    Ambmb0S  = (\[CapitalSigma]1/\[CapitalSigma]*vmb + 2*vmb1)*vmb + Cdmbmb;
    Ambmbt\[Phi]S = (St\[Phi]mb + exp1*vmb)*vmb;
    AmbmbrS  = (Srmb + rp1*vmb)*vmb;
    Ambmb\[Theta]S  = (S\[Theta]mb - zp1/sin\[Theta]p*vmb)*vmb;
      sumPlus0  += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RInrp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRInrp + Cmmbmb*fmbmb2*ddRInrp)*exp\[Theta]^{1,-1})]; (* Totral of all quadrants *)
      sumMinus0 += Total[\[CapitalSigma]*(((Cmnn*fnn0+Cmnmb*fnmb0+Cmmbmb*fmbmb0)*RUprp - (Cmnmb*fnmb1+Cmmbmb*fmbmb1)*dRUprp + Cmmbmb*fmbmb2*ddRUprp)*exp\[Theta]^{1,-1})]; 
      sumPlus1  += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RInrp*fnn0 + AnnrS*(dRInrp*fnn0 + RInrp*dfnn0dr) + Ann\[Theta]S*RInrp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RInrp*fnmb0 + AnmbrS*( RInrp*dfnmb0dr +  dRInrp*fnmb0) + Anmb\[Theta]S* RInrp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRInrp*fnmb1 + AnmbrS*(dRInrp*dfnmb1dr + ddRInrp*fnmb1) + Anmb\[Theta]S*dRInrp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RInrp*fmbmb0 + AmbmbrS*(  dRInrp*fmbmb0 +   RInrp*dfmbmb0dr) + Ambmb\[Theta]S*  RInrp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRInrp*fmbmb1 + AmbmbrS*( ddRInrp*fmbmb1 +  dRInrp*dfmbmb1dr) + Ambmb\[Theta]S* dRInrp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRInrp*fmbmb2 + AmbmbrS*(dddRInrp*fmbmb2 + ddRInrp*dfmbmb2dr) + Ambmb\[Theta]S*ddRInrp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RInrp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRInrp1 + Cmmbmb*fmbmb2*ddRInrp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RInrp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRInrp + Cmmbmb*fmbmb21*ddRInrp)*exp\[Theta]^{1,-1})]; 
      sumMinus1 += Total[\[CapitalSigma]*((
        ((Ann0S + Annt\[Phi]S)*RUprp*fnn0 + AnnrS*(dRUprp*fnn0 + RUprp*dfnn0dr) + Ann\[Theta]S*RUprp*dfnn0d\[Theta]) +
        ((Anmb0S + Anmbt\[Phi]S)* RUprp*fnmb0 + AnmbrS*( RUprp*dfnmb0dr +  dRUprp*fnmb0) + Anmb\[Theta]S* RUprp*dfnmb0d\[Theta]) - 
        ((Anmb0S + Anmbt\[Phi]S)*dRUprp*fnmb1 + AnmbrS*(dRUprp*dfnmb1dr + ddRUprp*fnmb1) + Anmb\[Theta]S*dRUprp*dfnmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*  RUprp*fmbmb0 + AmbmbrS*(  dRUprp*fmbmb0 +   RUprp*dfmbmb0dr) + Ambmb\[Theta]S*  RUprp*dfmbmb0d\[Theta]) -
        ((Ambmb0S + Ambmbt\[Phi]S)* dRUprp*fmbmb1 + AmbmbrS*( ddRUprp*fmbmb1 +  dRUprp*dfmbmb1dr) + Ambmb\[Theta]S* dRUprp*dfmbmb1d\[Theta]) +
        ((Ambmb0S + Ambmbt\[Phi]S)*ddRUprp*fmbmb2 + AmbmbrS*(dddRUprp*fmbmb2 + ddRUprp*dfmbmb2dr) + Ambmb\[Theta]S*ddRUprp*dfmbmb2d\[Theta]))*exp\[Theta]^{1,-1}) + 
        \[CapitalSigma]*(((Cmnn*fnn0 + Cmnmb*fnmb0 + Cmmbmb*fmbmb0)*RUprp1 - (Cmnmb*fnmb1 + Cmmbmb*fmbmb1)*dRUprp1 + Cmmbmb*fmbmb2*ddRUprp1
          + (Cmnn*fnn01 + Cmnmb*fnmb01 + Cmmbmb*fmbmb01)*RUprp - (Cmnmb*fnmb11 + Cmmbmb*fmbmb11)*dRUprp + Cmmbmb*fmbmb21*ddRUprp)*exp\[Theta]^{1,-1})]; 
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
    "steps\[Theta]"->steps\[Theta]
  |> (* l, m, k, n, \[Omega], C^+, C^-, \[Alpha], S(\[Pi]/2), dE^\[Infinity]/dt, dE^H/dt, Subscript[dJ, z]^\[Infinity]/dt, Subscript[dJ, z]^H/dt *)
]*)


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];

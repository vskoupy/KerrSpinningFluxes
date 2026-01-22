(* ::Package:: *)

(* ::Section::Closed:: *)
(*Begin package*)


BeginPackage["SpinningOrbit`",
  {
    "KerrGeodesics`KerrGeoOrbit`",
    "KerrGeodesics`OrbitalFrequencies`",
    "KerrGeodesics`ConstantsOfMotion`"
  }
];


KerrSpinOrbit::usage = "KerrSpinOrbit[a, p, e, x, spar] calculates trajectory of a spinning particle"
KerrSpinLinearParts::usage = "KerrSpinLinearParts[a, p, e, x] calculates linear parts of the constants, shifted constants, parameters, frequencies and actions in different gauges"
KerrSpinOrbitCorrection::usage = "KerrSpinOrbitCorrection[a, p, e, x, nmax, kmax] calculates linear correction to the trajectory";
KerrSpinOrbitNum::usage = "KerrSpinOrbitNum[KerrSpinOrbitCorrection, spar] calculates trajectory of spinning particle"
KerrSpinOrbitCorrectionSpherical::usage = "KerrSpinOrbitCorrectionSpherical[a, p, e, x, kmax] calculates linear correction to the spherical trajectory";
KerrSpinOrbitSpherical::usage = "KerrSpinOrbitSpherical[KerrSpinOrbitCorrection, spar] calculates nearly spherical trajectory of spinning particle"
KerrSpinOrbitEquatorial::usage = "KerrSpinOrbitEquatorial[orbit] returns eccentric equatorial orbit"


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Spinning trajectory*)


(* ::Subsection::Closed:: *)
(*Analytical generic orbit*)


(* ::Subsubsection::Closed:: *)
(*Linear parts*)


\[CapitalOmega]1CTildeFun[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]t,K},
\[CapitalUpsilon]\[Tau]=KerrGeodesics`OrbitalFrequencies`Private`KerrGeoProperFrequencyFactor[a,p,e,x];
\[CapitalUpsilon]t=KerrGeoFrequencies[a,p,e,x,"Time"->"Mino"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"];
K=KerrGeoCarterConstant[a,p,e,x]+(KerrGeoAngularMomentum[a,p,e,x]-a KerrGeoEnergy[a,p,e,x])^2;
Values[KerrGeoFrequencies[a,p,e,x] 3 \[CapitalUpsilon]\[Tau]/(2 \[CapitalUpsilon]t Sqrt[K])]
]
C1CTildeFun[a_,p_,e_,x_]:=Module[{CoM,En,Lz,Q,K},
CoM=KerrGeoConstantsOfMotion[a,p,e,x];
En=CoM["\[ScriptCapitalE]"];
Lz=CoM["\[ScriptCapitalL]"];
Q=CoM["\[ScriptCapitalQ]"];
K=Q+(Lz-a En)^2;
-{(1-En^2)/2,(a-Lz En/2),(3a (Lz-a En)-K En)}/Sqrt[K]
]
pex1CTildeFun[a_,p_,e_,x_]:=Module[{En,Lz,K,Q,r1,r2,r3,r4,z1,az2,kzsq,krsq,az1Int\[CapitalSigma],int\[CapitalSigma]r1,int\[CapitalSigma]r2,\[Delta]r1,\[Delta]r2,\[Delta]z1,\[Delta]p,\[Delta]e,\[Delta]x},
En=KerrGeoEnergy[a,p,e,x];
Lz=KerrGeoAngularMomentum[a,p,e,x,En];
Q=KerrGeoCarterConstant[a,p,e,x,En,Lz];
K=Q+(Lz-a En)^2;
{r1,r2,r3,r4}=KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En,Q];
z1=Sqrt[1-x^2];
az2=Sqrt[a^2 (1-z1^2)+(Lz^2+Q)/(1-En^2)];
kzsq=a^2 z1^2/az2^2;
krsq=(r1-r2) (r3-r4)/((r1-r3) (r2-r4));
int\[CapitalSigma]r1=EllipticPi[-a^2 z1^2/r1^2,kzsq]/EllipticK[kzsq]/r1^2;
int\[CapitalSigma]r2=EllipticPi[-a^2 z1^2/r2^2,kzsq]/EllipticK[kzsq]/r2^2;
az1Int\[CapitalSigma]=((a z1)/(r3^2+a^2 z1^2)+2 Re[((I (-r2+r3))/(2 (r2+I a z1) (r3+I a z1))) EllipticPi[((r1-r2) (r3+I a z1))/((r1-r3) (r2+I a z1)),krsq]]/EllipticK[krsq]);
\[Delta]r1=-r1 ((r1^2+a^2) En-a Lz) int\[CapitalSigma]r1/Sqrt[K];
\[Delta]r2=-r2 ((r2^2+a^2) En-a Lz) int\[CapitalSigma]r2/Sqrt[K];
\[Delta]z1=-(Lz-a En x^2) az1Int\[CapitalSigma]/Sqrt[K];
\[Delta]p=(2 r2^2 \[Delta]r1)/(r1+r2)^2+(2 r1^2 \[Delta]r2)/(r1+r2)^2;
\[Delta]e=(2 r2 \[Delta]r1)/(r1+r2)^2-(2 r1 \[Delta]r2)/(r1+r2)^2;
\[Delta]x=-((z1 \[Delta]z1)/x);
{\[Delta]p,\[Delta]e,\[Delta]x}
]
J1CFun[a_,p_,e_,x_]:=Module[{En,Lz,Q,K,az2,z1,Jz1,r1,r2,r3,r4,krsq,Jr1},
En=KerrGeoEnergy[a,p,e,x];
Lz=KerrGeoAngularMomentum[a,p,e,x,En];
Q=KerrGeoCarterConstant[a,p,e,x,En,Lz];
K=Q+(Lz-a En)^2;
{r1,r2,r3,r4}=KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En,Q];
krsq=((r1-r2) (r3-r4))/((r1-r3) (r2-r4));
Jr1=-2/\[Pi] Sqrt[K/((1-En^2)(r1-r3)(r2-r4))]((a^2 En - a Lz+En r3^2)/(K+r3^2) EllipticK[krsq]+((r2-r3) (-a^2 En+a Lz+En K))/Sqrt[K] Im[1/((r2-I Sqrt[K]) ( r3-I Sqrt[K])) EllipticPi[((r1-r2) (Sqrt[K]+I r3))/((Sqrt[K]+I r2) (r1-r3)),krsq]]);
z1=Sqrt[1-x^2];
az2=Sqrt[a^2 (1-z1^2)+(Lz^2+Q)/(1-En^2)];
Jz1=2/\[Pi] 1/(az2 Sqrt[K] Sqrt[1 - En^2]) (En K EllipticK[((a z1)/az2)^2] + (En(a^2 - K) - a Lz)EllipticPi[(a^2 z1^2)/K,((a z1)/az2)^2]);
{Jr1,Jz1,0}
]


(* ::Subsubsection::Closed:: *)
(*Jacobians*)


d\[CapitalOmega]JdCFun[a_,p_,e_,x_]:=Module[{M=1,En,Lz,Q,K,r1,r2,r3,r4,krsq,Kkr,Ekr,Pi\[Alpha]rkr,rplus,rminus,Pi\[Rho]pluskr,Pi\[Rho]minuskr,\[ScriptCapitalC],dJrdK,dJrdEn,dJrdLz,int1,int2,int3,int4,d2JrdK2,d2JrdKdEn,d2JrdKdLz,intr,
intr2,d2JrdEn2,d2JrdEndLz,d2JrdLz2,\[CapitalUpsilon]r,d\[CapitalUpsilon]rdK,d\[CapitalUpsilon]rdEn,d\[CapitalUpsilon]rdLz,\[CapitalUpsilon]tr,\[CapitalUpsilon]\[Phi]r,d\[CapitalUpsilon]trdEn,d\[CapitalUpsilon]trdLz,d\[CapitalUpsilon]trdK,d\[CapitalUpsilon]\[Phi]rdEn,d\[CapitalUpsilon]\[Phi]rdLz,d\[CapitalUpsilon]\[Phi]rdK,az2,z1,kzsq,Kkz,Ekz,Dkz,Pikz,dJzdK,dJzdEn,dJzdLz,d2JzdKdEn,d2JzdKdLz,d2JzdK2,d2JzdEn2,d2JzdEndLz,d2JzdLz2,\[CapitalUpsilon]z,d\[CapitalUpsilon]zdK,d\[CapitalUpsilon]zdEn,d\[CapitalUpsilon]zdLz,\[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]z,d\[CapitalUpsilon]tzdEn,d\[CapitalUpsilon]tzdLz,d\[CapitalUpsilon]tzdK,d\[CapitalUpsilon]\[Phi]zdEn,d\[CapitalUpsilon]\[Phi]zdLz,d\[CapitalUpsilon]\[Phi]zdK,\[CapitalUpsilon]t,\[CapitalUpsilon]\[Phi],d\[CapitalUpsilon]tdK,d\[CapitalUpsilon]tdEn,d\[CapitalUpsilon]tdLz,d\[CapitalUpsilon]\[Phi]dK,d\[CapitalUpsilon]\[Phi]dEn,d\[CapitalUpsilon]\[Phi]dLz,d\[CapitalOmega]rdK,d\[CapitalOmega]rdEn,d\[CapitalOmega]rdLz,d\[CapitalOmega]zdK,d\[CapitalOmega]zdEn,d\[CapitalOmega]zdLz,d\[CapitalOmega]\[Phi]dK,d\[CapitalOmega]\[Phi]dEn,d\[CapitalOmega]\[Phi]dLz,Iz2,Iz1},
En=KerrGeoEnergy[a,p,e,x];
Lz=KerrGeoAngularMomentum[a,p,e,x,En];
Q=KerrGeoCarterConstant[a,p,e,x,En,Lz];
K=Q+(Lz-a En)^2;
{r1,r2,r3,r4}=KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En,Q];
krsq=((r1-r2) (r3-r4))/((r1-r3) (r2-r4));
Kkr=EllipticK[krsq];
Ekr=EllipticE[krsq];
Pi\[Alpha]rkr=EllipticPi[(r1-r2)/(r1-r3),krsq];
rminus=1-Sqrt[1-a^2];
rplus=1+Sqrt[1-a^2];
Pi\[Rho]pluskr=EllipticPi[(r1-r2)/(r1-r3)*(r3-rplus)/(r2-rplus),krsq];
Pi\[Rho]minuskr=EllipticPi[(r1-r2)/(r1-r3)*(r3-rminus)/(r2-rminus),krsq];
\[ScriptCapitalC]=2/Sqrt[(1-En^2) (r1-r3) (r2-r4)];
dJrdK=-((\[ScriptCapitalC] Kkr)/(2\[Pi]));
dJrdEn=\[ScriptCapitalC]/\[Pi] (Kkr ((4+a^2) En-a Lz)+En (1/2 (Kkr (-r1 r2+r3 (4+r1+r2+r3))+Ekr (r1-r3) (r2-r4)+Pi\[Alpha]rkr (r2-r3) (4+r1+r2+r3+r4))+2/(-rminus+rplus) (-(((Kkr-(Pi\[Rho]minuskr (r2-r3))/(r2-rminus)) (-2 a^2+(4-(a Lz)/En) rminus))/(r3-rminus))+((Kkr-(Pi\[Rho]pluskr (r2-r3))/(r2-rplus)) (-2 a^2+(4-(a Lz)/En) rplus))/(r3-rplus))));
dJrdLz=-(\[ScriptCapitalC]/\[Pi]) (a En Kkr+a/(-rminus+rplus) (-(((Kkr-(Pi\[Rho]minuskr (r2-r3))/(r2-rminus)) (-a Lz+2 En rminus))/(r3-rminus))+((Kkr-(Pi\[Rho]pluskr (r2-r3))/(r2-rplus)) (-a Lz+2 En rplus))/(r3-rplus)));
int1=((r2-r4) Ekr+(-r1+r4) Kkr)/((r1-r2) (r1-r4));
int2=((-r1+r3) Ekr+(r2-r3) Kkr)/((r1-r2) (r2-r3));
int3=((r2-r4) Ekr+(-r2+r3) Kkr)/((r2-r3) (r3-r4));
int4=((r1-r3) Ekr+(-r1+r4) Kkr)/((r3-r4) (-r1+r4));
d2JrdK2=\[ScriptCapitalC]/(4Pi (1-En^2)) ((a^2-2 M r1+r1^2)/((r1-r2) (r1-r3) (r1-r4)) int1-(a^2-2 M r2+r2^2)/((r1-r2) (r2-r3) (r2-r4)) int2-(a^2-2 M r3+r3^2)/((r1-r3) (-r2+r3) (r3-r4)) int3-(a^2-2 M r4+r4^2)/((r1-r4) (-r2+r4) (-r3+r4)) int4);
d2JrdKdEn=-\[ScriptCapitalC]/(2Pi (1-En^2)) (En Kkr+((a^2+r1^2) (a^2 En-a Lz+En r1^2))/((r1-r2) (r1-r3) (r1-r4))  int1-((a^2+r2^2) (a^2 En-a Lz+En r2^2))/((r1-r2) (r2-r3) (r2-r4)) int2-((a^2+r3^2) (a^2 En-a Lz+En r3^2))/((r1-r3) (-r2+r3) (r3-r4)) int3-((a^2+r4^2) (a^2 En-a Lz+En r4^2))/((r1-r4) (-r2+r4) (-r3+r4)) int4);
d2JrdKdLz=a \[ScriptCapitalC]/(2Pi (1-En^2)) ((a^2 En-a Lz+En r1^2)/((r1-r2) (r1-r3) (r1-r4)) int1-(a^2 En-a Lz+En r2^2)/((r1-r2) (r2-r3) (r2-r4)) int2-(a^2 En-a Lz+En r3^2)/((r1-r3) (-r2+r3) (r3-r4)) int3-(a^2 En-a Lz+En r4^2)/((r1-r4) (-r2+r4) (-r3+r4)) int4);
intr=r3 Kkr+(r2-r3) Pi\[Alpha]rkr;
intr2=1/2 (r1 (-r2+r3)+r3 (r2+r3)) Kkr+1/2 (r1-r3) (r2-r4) Ekr+1/2 (r2-r3) (r1+r2+r3+r4) Pi\[Alpha]rkr;
d2JrdEn2=\[ScriptCapitalC]/(Pi (1-En^2)) ((4+(-1+En^2) (-a^2+En^2 K+2 a En Lz))/(-1+En^2)^2 Kkr+2/(1-En^2) intr+intr2+((a^2+r1^2)^2 (K+r1^2))/((r1-r2) (r1-r3) (r1-r4)) int1-((a^2+r2^2)^2 (K+r2^2))/((r1-r2) (r2-r3) (r2-r4)) int2-((a^2+r3^2)^2 (K+r3^2))/((r1-r3) (-r2+r3) (r3-r4)) int3-((a^2+r4^2)^2 (K+r4^2))/((r1-r4) (-r2+r4) (-r3+r4)) int4);
d2JrdEndLz=\[ScriptCapitalC] a/(Pi (1-En^2)) (-Kkr-(int1 (a^2+r1^2) (K+r1^2))/((r1-r2) (r1-r3) (r1-r4))+(int2 (a^2+r2^2) (K+r2^2))/((r1-r2) (r2-r3) (r2-r4))+(int3 (a^2+r3^2) (K+r3^2))/((r1-r3) (-r2+r3) (r3-r4))+(int4 (a^2+r4^2) (K+r4^2))/((r1-r4) (-r2+r4) (-r3+r4)));
d2JrdLz2=\[ScriptCapitalC] a^2/(Pi (1-En^2)) ((K+r1^2)/((r1-r2) (r1-r3) (r1-r4)) int1+(-K-r2^2)/((r1-r2) (r2-r3) (r2-r4)) int2+(-K-r3^2)/((r1-r3) (-r2+r3) (r3-r4)) int3+(-K-r4^2)/((r1-r4) (-r2+r4) (-r3+r4)) int4);
\[CapitalUpsilon]r=-1/(2 dJrdK);
d\[CapitalUpsilon]rdK=d2JrdK2/(2 dJrdK^2);
d\[CapitalUpsilon]rdEn=d2JrdKdEn/(2 dJrdK^2);
d\[CapitalUpsilon]rdLz=d2JrdKdLz/(2 dJrdK^2);
\[CapitalUpsilon]tr=-dJrdEn/(2 dJrdK);
\[CapitalUpsilon]\[Phi]r=dJrdLz/(2 dJrdK);
d\[CapitalUpsilon]trdK=-d2JrdKdEn/(2 dJrdK)+dJrdEn/(2 dJrdK^2) d2JrdK2;
d\[CapitalUpsilon]trdEn=-d2JrdEn2/(2 dJrdK)+dJrdEn/(2 dJrdK^2) d2JrdKdEn;
d\[CapitalUpsilon]trdLz=-d2JrdEndLz/(2 dJrdK)+dJrdEn/(2 dJrdK^2) d2JrdKdLz;
d\[CapitalUpsilon]\[Phi]rdK=d2JrdKdLz/(2 dJrdK)-dJrdLz/(2 dJrdK^2) d2JrdK2;
d\[CapitalUpsilon]\[Phi]rdEn=d2JrdEndLz/(2 dJrdK)-dJrdLz/(2 dJrdK^2) d2JrdKdEn;
d\[CapitalUpsilon]\[Phi]rdLz=d2JrdLz2/(2 dJrdK)-dJrdLz/(2 dJrdK^2) d2JrdKdLz;
(********************************************************)
z1=Sqrt[1-x^2];
az2=Sqrt[a^2 (1-z1^2)+(Lz^2+Q)/(1-En^2)];
kzsq=(a z1/az2)^2;
Kkz=EllipticK[kzsq];
Ekz=EllipticE[kzsq];
Dkz=1/3CarlsonRD[0,1-kzsq,1];
Pikz=EllipticPi[z1^2,kzsq];
dJzdK=Kkz/(Sqrt[1-En^2] \[Pi] az2);
dJzdEn=1/(Sqrt[1-En^2] \[Pi] az2) 2 (-En az2^2 Ekz+(a Lz-En (a^2-az2^2)) Kkz);
dJzdLz=-((2 (-a En Kkz+Lz Pikz))/(Sqrt[1-En^2] \[Pi] az2));
Iz2=(-(Ekz/(1-kzsq)));
Iz1=((-Dkz+Kkz)/(-1+kzsq));
d2JzdK2=-1/(2Pi az2^3 (1-En^2)^(3/2) (az2^2-a^2 z1^2)) ((a^2-az2^2)Iz2+(1-z1^2)a^2 Iz1);
d2JzdKdEn=1/(Pi az2 (1-En^2)^(3/2)) (En Kkz-1/((az2^2-a^2 z1^2) az2^2) ((a^2-az2^2) (a Lz+En (-a^2+az2^2))Iz2+a^3 (1-z1^2) (Lz+a En (-1+z1^2)) Iz1));
d2JzdKdLz=1/(Pi az2^3 (1-En^2)^(3/2) (az2^2-a^2 z1^2)) ((a^2 Lz+a En (-a^2+az2^2))Iz2+(Lz+a En (-1+z1^2))a^2 Iz1);
d2JzdEn2=2/(Pi az2 (1-En^2)^(3/2)) (-(K- (-2 a^2+az2^2+a^2 z1^2)) Kkz+az2^2 (Kkz-Ekz)-1/((az2^2-a^2 z1^2) az2^2) ((a^2-az2^2)^2 (K-az2^2)Iz2+a^4 (1-z1^2)^2 (K-a^2 z1^2) Iz1));
d2JzdEndLz=2 a/(Pi az2 (1-En^2)^(3/2)) (Kkz+1/((az2^2-a^2 z1^2) az2^2) ((a^2-az2^2) (K-az2^2)Iz2+a^2 (1-z1^2) (K-a^2 z1^2) Iz1));
d2JzdLz2=-2 a^2/(Pi az2^3 (1-En^2)^(3/2) (az2^2-a^2 z1^2)) ((K-az2^2)Iz2+(K-a^2 z1^2) Iz1);
\[CapitalUpsilon]z=1/(2 dJzdK);
d\[CapitalUpsilon]zdK=-d2JzdK2/(2 dJzdK^2);
d\[CapitalUpsilon]zdEn=-d2JzdKdEn/(2 dJzdK^2);
d\[CapitalUpsilon]zdLz=-d2JzdKdLz/(2 dJzdK^2);
\[CapitalUpsilon]tz=dJzdEn/(2 dJzdK);
\[CapitalUpsilon]\[Phi]z=-dJzdLz/(2 dJzdK);
d\[CapitalUpsilon]tzdK=d2JzdKdEn/(2 dJzdK)-dJzdEn/(2 dJzdK^2) d2JzdK2;
d\[CapitalUpsilon]tzdEn=d2JzdEn2/(2 dJzdK)-dJzdEn/(2 dJzdK^2) d2JzdKdEn;
d\[CapitalUpsilon]tzdLz=d2JzdEndLz/(2 dJzdK)-dJzdEn/(2 dJzdK^2) d2JzdKdLz;
d\[CapitalUpsilon]\[Phi]zdK=-d2JzdKdLz/(2 dJzdK)+dJzdLz/(2 dJzdK^2) d2JzdK2;
d\[CapitalUpsilon]\[Phi]zdEn=-d2JzdEndLz/(2 dJzdK)+dJzdLz/(2 dJzdK^2) d2JzdKdEn;
d\[CapitalUpsilon]\[Phi]zdLz=-d2JzdLz2/(2 dJzdK)+dJzdLz/(2 dJzdK^2) d2JzdKdLz;
(********************************************************)
\[CapitalUpsilon]t=\[CapitalUpsilon]tr+\[CapitalUpsilon]tz;
\[CapitalUpsilon]\[Phi]=\[CapitalUpsilon]\[Phi]r+\[CapitalUpsilon]\[Phi]z;
d\[CapitalUpsilon]tdK=d\[CapitalUpsilon]trdK+d\[CapitalUpsilon]tzdK;
d\[CapitalUpsilon]tdEn=d\[CapitalUpsilon]trdEn+d\[CapitalUpsilon]tzdEn;
d\[CapitalUpsilon]tdLz=d\[CapitalUpsilon]trdLz+d\[CapitalUpsilon]tzdLz;
d\[CapitalUpsilon]\[Phi]dK=d\[CapitalUpsilon]\[Phi]rdK+d\[CapitalUpsilon]\[Phi]zdK;
d\[CapitalUpsilon]\[Phi]dEn=d\[CapitalUpsilon]\[Phi]rdEn+d\[CapitalUpsilon]\[Phi]zdEn;
d\[CapitalUpsilon]\[Phi]dLz=d\[CapitalUpsilon]\[Phi]rdLz+d\[CapitalUpsilon]\[Phi]zdLz;
d\[CapitalOmega]rdK=d\[CapitalUpsilon]rdK/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdK \[CapitalUpsilon]r/\[CapitalUpsilon]t^2;
d\[CapitalOmega]rdEn=d\[CapitalUpsilon]rdEn/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdEn \[CapitalUpsilon]r/\[CapitalUpsilon]t^2;
d\[CapitalOmega]rdLz=d\[CapitalUpsilon]rdLz/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdLz \[CapitalUpsilon]r/\[CapitalUpsilon]t^2;
d\[CapitalOmega]zdK=d\[CapitalUpsilon]zdK/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdK \[CapitalUpsilon]z/\[CapitalUpsilon]t^2;
d\[CapitalOmega]zdEn=d\[CapitalUpsilon]zdEn/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdEn \[CapitalUpsilon]z/\[CapitalUpsilon]t^2;
d\[CapitalOmega]zdLz=d\[CapitalUpsilon]zdLz/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdLz \[CapitalUpsilon]z/\[CapitalUpsilon]t^2;
d\[CapitalOmega]\[Phi]dK=d\[CapitalUpsilon]\[Phi]dK/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdK \[CapitalUpsilon]\[Phi]/\[CapitalUpsilon]t^2;
d\[CapitalOmega]\[Phi]dEn=d\[CapitalUpsilon]\[Phi]dEn/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdEn \[CapitalUpsilon]\[Phi]/\[CapitalUpsilon]t^2;
d\[CapitalOmega]\[Phi]dLz=d\[CapitalUpsilon]\[Phi]dLz/\[CapitalUpsilon]t-d\[CapitalUpsilon]tdLz \[CapitalUpsilon]\[Phi]/\[CapitalUpsilon]t^2;
{{{d\[CapitalOmega]rdEn,d\[CapitalOmega]rdLz,d\[CapitalOmega]rdK},{d\[CapitalOmega]zdEn,d\[CapitalOmega]zdLz,d\[CapitalOmega]zdK},{d\[CapitalOmega]\[Phi]dEn,d\[CapitalOmega]\[Phi]dLz,d\[CapitalOmega]\[Phi]dK}},
{{dJrdEn,dJrdLz,dJrdK},{dJzdEn,dJzdLz,dJzdK},{0,1,0}}}
]
dr12z1sqdCFun[a_,p_,e_,x_]:=Module[{dRdr,dRdEn,dRdLz,dRdK,dZdzsq,dZdEn,dZdLz,dZdK,CoM,En,Lz,Q,K,r1,r2,z1,\[CapitalDelta]},
\[CapitalDelta][r_]:=r^2-2r+a^2;
dRdr[r_]:=4 En r (-a Lz+En (a^2+r^2))-2 r \[CapitalDelta][r]-((-a En+Lz)^2+Q+r^2) Derivative[1][\[CapitalDelta]][r];
dRdEn[r_]:=2 (a^2+r^2) (a^2 En-a Lz+En r^2);
dRdLz[r_]:=2 a^2 Lz-2 a En (a^2+r^2);
dRdK[r_]:=-\[CapitalDelta][r];
dZdzsq[z_]:=-(Lz^2+Q+a^2 (-1+En^2) (-1+2 z^2));
dZdEn[z_]:=2 a (1-z^2) (Lz-a En (1-z^2));
dZdLz[z_]:=-2 (Lz-a En (1-z^2));
dZdK[z_]:=1-z^2;
r1=p/(1-e);
r2=p/(1+e);
z1=Sqrt[1-x^2];
CoM=KerrGeoConstantsOfMotion[a,p,e,x];
En=CoM["\[ScriptCapitalE]"];
Lz=CoM["\[ScriptCapitalL]"];
Q=CoM["\[ScriptCapitalQ]"];
K=Q+(Lz-a En)^2;
-{{dRdEn[r1],dRdLz[r1],dRdK[r1]}/dRdr[r1],{dRdEn[r2],dRdLz[r2],dRdK[r2]}/dRdr[r2],{dZdEn[z1],dZdLz[z1],dZdK[z1]}/dZdzsq[z1]}
]
dpexdr12z1sqFun[a_,p_,e_,x_]:={{1/2 (-1+e)^2,1/2 (1+e)^2,0},{((-1+e)^2 (1+e))/(2 p),((-1+e) (1+e)^2)/(2 p),0},{0,0,-(1/(2 x))}}
dpexdCFun[a_,p_,e_,x_]:=dpexdr12z1sqFun[a,p,e,x] . dr12z1sqdCFun[a,p,e,x]


(* ::Subsubsection::Closed:: *)
(*Master function for linear parts*)


KerrSpinLinearParts[a_,p_,e_,x_]:=Module[{\[CapitalOmega]1CTilde,C1CTilde,pex1CTilde,J1C,d\[CapitalOmega]dC,dJdC,dpexdC,CTilde1C,CTilde1pex,CTilde1\[CapitalOmega],\[CapitalOmega]1C,\[CapitalOmega]1pex,C1pex,pex1C,C1\[CapitalOmega],pex1\[CapitalOmega],J1CTilde,J1pex,J1\[CapitalOmega],\[CapitalOmega]1J,C1J,CTilde1J,pex1J},
\[CapitalOmega]1CTilde=\[CapitalOmega]1CTildeFun[a,p,e,x];
C1CTilde=C1CTildeFun[a,p,e,x];
pex1CTilde=pex1CTildeFun[a,p,e,x];
J1C=J1CFun[a,p,e,x];
{d\[CapitalOmega]dC,dJdC}=d\[CapitalOmega]JdCFun[a,p,e,x];
dpexdC=dpexdCFun[a,p,e,x];
CTilde1C=-C1CTilde;
CTilde1pex=-Inverse[dpexdC] . pex1CTilde;
CTilde1\[CapitalOmega]=-Inverse[d\[CapitalOmega]dC] . \[CapitalOmega]1CTilde;
\[CapitalOmega]1C=\[CapitalOmega]1CTilde-d\[CapitalOmega]dC . C1CTilde;
\[CapitalOmega]1pex=\[CapitalOmega]1CTilde-d\[CapitalOmega]dC . Inverse[dpexdC] . pex1CTilde;
C1pex=C1CTilde-Inverse[dpexdC] . pex1CTilde;
pex1C=pex1CTilde-dpexdC . C1CTilde;
C1\[CapitalOmega]=C1CTilde-Inverse[d\[CapitalOmega]dC] . \[CapitalOmega]1CTilde;
pex1\[CapitalOmega]=pex1CTilde-dpexdC . Inverse[d\[CapitalOmega]dC] . \[CapitalOmega]1CTilde;
J1CTilde=J1C+dJdC . C1CTilde;
J1pex=J1C+dJdC . (C1CTilde-Inverse[dpexdC] . pex1CTilde);
J1\[CapitalOmega]=J1C+dJdC . (C1CTilde-Inverse[d\[CapitalOmega]dC] . \[CapitalOmega]1CTilde);
\[CapitalOmega]1J=\[CapitalOmega]1CTilde-d\[CapitalOmega]dC . (C1CTilde+Inverse[dJdC] . J1C);
C1J=-Inverse[dJdC] . J1C;
CTilde1J=-C1CTilde-Inverse[dJdC] . J1C;
pex1J=pex1CTilde-dpexdC . (C1CTilde+Inverse[dJdC] . J1C);
<|
"\[CapitalOmega]1CTilde"->\[CapitalOmega]1CTilde,
"C1CTilde"->C1CTilde,
"pex1CTilde"->pex1CTilde,
"J1C"->J1C,
"CTilde1C"->CTilde1C,
"CTilde1pex"->CTilde1pex,
"CTilde1\[CapitalOmega]"->CTilde1\[CapitalOmega],
"\[CapitalOmega]1C"->\[CapitalOmega]1C,
"\[CapitalOmega]1pex"->\[CapitalOmega]1pex,
"C1pex"->C1pex,
"pex1C"->pex1C,
"C1\[CapitalOmega]"->C1\[CapitalOmega],
"pex1\[CapitalOmega]"->pex1\[CapitalOmega],
"J1CTilde"->J1CTilde,
"J1pex"->J1pex,
"J1\[CapitalOmega]"->J1\[CapitalOmega],
"\[CapitalOmega]1J"->\[CapitalOmega]1J,
"C1J"->C1J,
"CTilde1J"->CTilde1J,
"pex1J"->pex1J,
"d\[CapitalOmega]dC"->d\[CapitalOmega]dC,
"dJdC"->dJdC,
"dpexdC"->dpexdC
|>
]


(* ::Subsubsection::Closed:: *)
(*Auxiliary functions*)


\[Delta]t3[r_, z_, Ur_, Uz_, K_, a_]:=(r (a^2+r^2) Ur - a^2 (a^2 + r (-2+r)) z Uz)/(Sqrt[K] (a^2 + r (-2+r)) (r^2+a^2 z^2));
\[Delta]r3[r_, z_, En_, Jz_, K_, a_]:=(r (a^2 En - a Jz+En r^2))/(Sqrt[K] (r^2+a^2 z^2));
\[Delta]z3[r_, z_, En_, Jz_, K_, a_]:=(a z (Jz - a En (1-z^2)))/(Sqrt[K] (r^2+a^2 z^2));
\[Delta]\[Phi]3[r_, z_, Ur_, Uz_, K_, a_]:=(a (r (1-z^2) Ur - (a^2 + r (-2+r)) z Uz))/(Sqrt[K] (a^2 + r (-2+r)) (1-z^2) (r^2+a^2 z^2));


KerrGeoProperTime[orbitGeo_]:=Module[{M=1,a,En,Lz,K,r1,r2,r3,r4,x,z1,z2,krsq,kzsq,
	Kkr,Kkz,\[CapitalUpsilon]r,\[CapitalUpsilon]z,rminus,rplus,hr,hplus,hminus,\[CapitalUpsilon]tr,\[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]r,\[CapitalUpsilon]\[Phi]z,\[CapitalUpsilon]\[Tau]r,\[CapitalUpsilon]\[Tau]z,\[CapitalUpsilon]\[Psi]r,\[CapitalUpsilon]\[Psi]z,\[CapitalUpsilon]t,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]\[Tau],\[CapitalUpsilon]\[Psi],rTilde,zTilde,tTilde,\[Phi]Tilde,\[Tau],\[Psi],thatr,thatz,\[Phi]hatr,\[Phi]hatz,\[Tau]hatr,\[Tau]hatz,\[Psi]hatr,\[Psi]hatz,\[CapitalDelta]trTilde,\[CapitalDelta]tzTilde,\[CapitalDelta]\[Phi]rTilde,\[CapitalDelta]\[Phi]zTilde,\[CapitalDelta]\[Tau]r,\[CapitalDelta]\[Tau]z,\[CapitalDelta]\[Psi]r,\[CapitalDelta]\[Psi]z,\[Delta]rpar,\[Delta]zpar,\[Delta]tpar,\[Delta]\[Phi]par,t,r,z,\[Phi]},
	a = orbitGeo["a"];
	En = orbitGeo["Energy"];
	Lz = orbitGeo["AngularMomentum"];
	K = orbitGeo["CarterConstant"]+(Lz-a*En)^2;
	{r1,r2,r3,r4}=orbitGeo["RadialRoots"];
	x = orbitGeo["Inclination"];
	z1=Sqrt[1-x^2];
	z2=Sqrt[(K-(Lz-a*En)^2)/(a^2*(1-En^2)*z1^2)];
	krsq=(r1-r2)*(r3-r4)/((r1-r3)*(r2-r4));
	kzsq=z1^2/z2^2;
	Kkr=EllipticK[krsq];
	Kkz=EllipticK[kzsq];
	\[CapitalUpsilon]r=Pi*Sqrt[(1-En^2)*(r1-r3)*(r2-r4)]/(2*Kkr);
	\[CapitalUpsilon]z=Pi*Sqrt[a^2*(1-En^2)]*z2/(2*Kkz);
	rminus=M-Sqrt[M^2-a^2];
	rplus=M+Sqrt[M^2-a^2];
	hr=(r1-r2)/(r1-r3);
	hplus=hr*(r3-rplus)/(r2-rplus);
	hminus=hr*(r3-rminus)/(r2-rminus);
	\[CapitalUpsilon]\[Tau]r=1/2*((r1 (-r2+r3)+r3 (r2+r3))+(r2-r3) (r1+r2+r3+r4)*EllipticPi[hr,krsq]/Kkr+(r1-r3) (r2-r4)*EllipticE[krsq]/Kkr);
	\[CapitalUpsilon]\[Tau]z=a^2*z2^2*(1-EllipticE[kzsq]/Kkz);
	\[CapitalUpsilon]\[Tau]=\[CapitalUpsilon]\[Tau]r+\[CapitalUpsilon]\[Tau]z;
	\[Tau]hatz[\[Xi]z_]:=-z2*a/Sqrt[1-En^2]*EllipticE[\[Xi]z,kzsq];
	\[Tau]hatr[\[Xi]r_]:=((r2-r3) (r1+r2+r3+r4)*EllipticPi[hr,\[Xi]r,krsq]+(r1-r3) (r2-r4)*EllipticE[\[Xi]r,krsq]-(r1-r2) (r2-r4)*Sin[\[Xi]r]*Cos[\[Xi]r]*Sqrt[1-krsq*Sin[\[Xi]r]^2]/(1-hr*Sin[\[Xi]r]^2))/Sqrt[(1-En^2)*(r1-r3)*(r2-r4)];
	\[CapitalDelta]\[Tau]r[qr_]:=\[Tau]hatr[JacobiAmplitude[Kkr*qr/Pi,krsq]]-\[Tau]hatr[Pi]/(2Pi)*qr;
	\[CapitalDelta]\[Tau]z[qz_]:=\[Tau]hatz[JacobiAmplitude[2Kkz*(qz+Pi/2)/Pi,kzsq]]-\[Tau]hatz[Pi]/Pi*(qz+Pi/2);
	\[Tau][q\[Tau]_,qr_,qz_]:=q\[Tau]+\[CapitalDelta]\[Tau]r[qr]+\[CapitalDelta]\[Tau]z[qz];
	<|
		"ProperTimeFrequency"->\[CapitalUpsilon]\[Tau],
		"ProperTime"->Function[{q\[Tau],qr,qz},\[Tau][q\[Tau],qr,qz]],
		"\[CapitalUpsilon]\[Tau]z" -> \[CapitalUpsilon]\[Tau]z
	|>
]


(* ::Subsubsection::Closed:: *)
(*KerrSpinOrbit function*)


Options[KerrSpinOrbit] = {"Gauge"->"CTilde"};
KerrSpinOrbit[a_, p_, e_, x_, spar_, OptionsPattern[]]:=Module[{linearParts,pTilde,eTilde,xTilde,orbitTilde,EnTilde,JzTilde,KTilde,
	sqrtKTilde,En,Jz,K,\[Tau]Tilde,frequencies,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalDelta]tr,\[CapitalDelta]tz,rTilde,zTilde,
	UrTilde,UzTilde,trajectoryTilde,\[Delta]x,\[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]z},
	linearParts = KerrSpinLinearParts[a,p,e,x];
	Switch[OptionValue["Gauge"],
		"CTilde",
		{pTilde,eTilde,xTilde} = {p,e,x};,
		"pex",
		{pTilde,eTilde,xTilde} = {p,e,x} - spar*linearParts["pex1CTilde"];,
		"C",
		{pTilde,eTilde,xTilde} = {p,e,x} - spar*linearParts["dpexdC"] . linearParts["C1CTilde"];,
		_,
		Print["Unknown gauge"];Return[$Failed];
	];
	orbitTilde = KerrGeodesics`KerrGeoOrbit`KerrGeoOrbit[a, pTilde, eTilde, xTilde, "Parametrization" -> "Phases", "Method" -> "Analytic"];
	EnTilde = orbitTilde["ConstantsOfMotion"]["\[ScriptCapitalE]"];
	JzTilde = orbitTilde["ConstantsOfMotion"]["\[ScriptCapitalL]"];
	KTilde = orbitTilde["ConstantsOfMotion"]["\[ScriptCapitalQ]"] + (JzTilde - a*EnTilde)^2;
	sqrtKTilde = Sqrt[KTilde];
	En = EnTilde - spar*(1 - EnTilde^2)/(2*sqrtKTilde);
	Jz = JzTilde - spar*(a - EnTilde*JzTilde/2)/(sqrtKTilde);
	K  = KTilde  - spar*(3*a*(JzTilde - a*EnTilde) - EnTilde*KTilde)/(sqrtKTilde);
	\[Tau]Tilde = KerrGeoProperTime[orbitTilde];
	frequencies = orbitTilde["Frequencies"];
	\[CapitalUpsilon]r = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"];
	\[CapitalUpsilon]z = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"];
	\[CapitalUpsilon]\[Phi] = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"];
	\[CapitalUpsilon]t = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"] - 3*spar/(2*sqrtKTilde)*\[Tau]Tilde["ProperTimeFrequency"];
	\[CapitalUpsilon]tz = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencyt\[Theta][a,pTilde,eTilde,xTilde,{EnTilde,JzTilde,KTilde-(JzTilde-a*EnTilde)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,pTilde,eTilde,xTilde]] + a*JzTilde - 3*spar/(2*sqrtKTilde)*\[Tau]Tilde["\[CapitalUpsilon]\[Tau]z"];
	\[CapitalUpsilon]\[Phi]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Phi]\[Theta][a,pTilde,eTilde,xTilde,{EnTilde,JzTilde,KTilde-(JzTilde-a*EnTilde)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,pTilde,eTilde,xTilde]] - a*EnTilde;
	\[CapitalDelta]tr = Function[qr, orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]tr"][qr] - 3*spar/(2*sqrtKTilde)*\[Tau]Tilde["ProperTime"][0,qr,0]];
	\[CapitalDelta]tz = Function[qz, orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][qz] - 3*spar/(2*sqrtKTilde)*\[Tau]Tilde["ProperTime"][0,0,qz]];
	rTilde = orbitTilde["Trajectory"][[2]];
	zTilde[qz_] := Cos[orbitTilde["Trajectory"][[3]][qz]];
	UrTilde[qr_]:=If[Mod[qr,2Pi]<Pi,1,-1]*Sqrt[((rTilde[qr]^2+a^2)*EnTilde - a*JzTilde)^2 - (rTilde[qr]^2 - 2*rTilde[qr] + a^2)*(KTilde + rTilde[qr]^2)];
	UzTilde[qz_]:=If[Mod[qz,2Pi]<Pi,-1,1]*Sqrt[-((1-zTilde[qz]^2)*a*EnTilde - JzTilde)^2 + (1-zTilde[qz]^2)*(KTilde - a^2*zTilde[qz]^2)];
	trajectoryTilde = {Function[{qt,qr,qz},orbitTilde["Trajectory"][[1]][qt,qr,qz] - 3*spar/(2*sqrtKTilde)*\[Tau]Tilde["ProperTime"][0,qr,qz]],
	rTilde,zTilde,orbitTilde["Trajectory"][[4]]};
	\[Delta]x = Function[{qr,qz},{
		\[Delta]t3[rTilde[qr], zTilde[qz], UrTilde[qr], UzTilde[qz], KTilde, a],
		\[Delta]r3[rTilde[qr], zTilde[qz], EnTilde, JzTilde, KTilde, a],
		\[Delta]z3[rTilde[qr], zTilde[qz], EnTilde, JzTilde, KTilde, a],
		\[Delta]\[Phi]3[rTilde[qr], zTilde[qz], UrTilde[qr], UzTilde[qz], KTilde, a]}];
	<|
		"a"->a,
		"p"->p,
		"e"->e,
		"x"->x,
		"spar"->spar,
		"Energy"->En,
		"AngularMomentum"->Jz,
		"CarterConstantK"->K,
		"EnTilde"->EnTilde,
		"JzTilde"->JzTilde,
		"KTilde"->KTilde,
		"pTilde"->pTilde,
		"eTilde"->eTilde,
		"xTilde"->xTilde,
		"Frequencies"-><|
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"->\[CapitalUpsilon]r,
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(z\)]\)"->\[CapitalUpsilon]z,
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"->\[CapitalUpsilon]\[Phi],
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"->\[CapitalUpsilon]t
		|>,
		"ProperTimeFrequency"->\[Tau]Tilde["ProperTimeFrequency"],
		"BLFrequencies"->{\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t,
		"TrajectoryDeltas"-><|
			"\[CapitalDelta]tr"->\[CapitalDelta]tr,
			"\[CapitalDelta]tz"->\[CapitalDelta]tz,
			"\[CapitalDelta]\[Phi]r"->orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"],
			"\[CapitalDelta]\[Phi]z"->orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]
		|>,
		"TrajectoryTilde"->trajectoryTilde,
		"\[Delta]x"->\[Delta]x,
		"FourVelocitiesTilde"->{Function[qr,UrTilde[qr]],Function[qz,UzTilde[qz]]},
		"\[CapitalUpsilon]tz"->\[CapitalUpsilon]tz,
		"\[CapitalUpsilon]\[Phi]z"->\[CapitalUpsilon]\[Phi]z
	|>
]


(*dr12z1dCOM[a_,p_,e_,x_,constants_]:=Module[{dRdr,dRdEn,dRdLz,dRdK,dZdz,dZdEn,dZdLz,dZdK,En,Lz,K,r1,r2,z1,\[CapitalDelta]},
\[CapitalDelta][r_]:=r^2-2r+a^2;
dRdr[r_]:=4 En r (-a Lz+En (a^2+r^2))-2 r \[CapitalDelta][r]-(K+r^2) Derivative[1][\[CapitalDelta]][r];
dRdEn[r_]:=2 (a^2+r^2) (a^2 En-a Lz+En r^2);
dRdLz[r_]:=2 a^2 Lz-2 a En (a^2+r^2);
dRdK[r_]:=-\[CapitalDelta][r];
dZdz[z_]:=2*z*(-(Lz^2+K-(Lz-a*En)^2+a^2 (-1+En^2) (-1+2 z^2)));
dZdEn[z_]:=2 a*(1-z^2)*(Lz-a*En*(1-z^2));
dZdLz[z_]:=-2*(Lz-a*En*(1-z^2));
dZdK[z_]:=1-z^2;
r1=p/(1-e);
r2=p/(1+e);
z1=Sqrt[1-x^2];
En=constants["\[ScriptCapitalE]"];
Lz=constants["\[ScriptCapitalL]"];
K=constants["\[ScriptCapitalQ]"]+(Lz-a*En)^2;
-{{dRdEn[r1]/dRdr[r1],dRdLz[r1]/dRdr[r1],dRdK[r1]/dRdr[r1]},{dRdEn[r2]/dRdr[r2],dRdLz[r2]/dRdr[r2],dRdK[r2]/dRdr[r2]},{dZdEn[z1]/dZdz[z1],dZdLz[z1]/dZdz[z1],dZdK[z1]/dZdz[z1]}}
]*)


(*KerrSpinOrbit[a_, p_, e_, x_, spar_]:=Module[{constants,En,Jz,K,sqrtK,
	\[Delta]En,\[Delta]Jz,\[Delta]K,EnTilde,JzTilde,KTilde,dr1dEn,dr1dLz,dr1dK,dr2dEn,dr2dLz,
	dr2dK,dz1dEn,dz1dLz,dz1dK,r1,r2,z1,\[Delta]r1,\[Delta]r2,\[Delta]z1,xTilde,r,pTilde,eTilde,
	orbitTilde,\[Tau]Tilde,frequencies,\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,\[CapitalDelta]tr,\[CapitalDelta]tz,rTilde,zTilde,
	UrTilde,UzTilde,Trajectoryg,\[Delta]x,\[CapitalUpsilon]tz,\[CapitalUpsilon]\[Phi]z},
	constants = KerrGeodesics`ConstantsOfMotion`KerrGeoConstantsOfMotion[a,p,e,x];
	En = constants["\[ScriptCapitalE]"];
	Jz = constants["\[ScriptCapitalL]"];
	K = constants["\[ScriptCapitalQ]"] + (Jz - a*En)^2;
	sqrtK = Sqrt[K];
	\[Delta]En = spar*(1 - En^2)/(2*sqrtK);
	\[Delta]Jz = spar*(a - En*Jz/2)/(sqrtK);
	\[Delta]K = spar*(3*a*(Jz - a*En) - En*K)/(sqrtK);
	EnTilde = En + \[Delta]En;
	JzTilde = Jz + \[Delta]Jz;
	KTilde = K + \[Delta]K;
	{{dr1dEn,dr1dLz,dr1dK},{dr2dEn,dr2dLz,dr2dK},{dz1dEn,dz1dLz,dz1dK}}=dr12z1dCOM[a,p,e,x,constants];
	r1 = p/(1-e);
	r2 = p/(1+e);
	z1 = Sqrt[1-x^2];
	{\[Delta]r1,\[Delta]r2,\[Delta]z1}={{dr1dEn,dr1dLz,dr1dK},{dr2dEn,dr2dLz,dr2dK},{dz1dEn,dz1dLz,dz1dK}} . {\[Delta]En,\[Delta]Jz,\[Delta]K};
	pTilde=2*(r1+\[Delta]r1)*(r2+\[Delta]r2)/(r1+\[Delta]r1+r2+\[Delta]r2);
	eTilde=((r1+\[Delta]r1)-(r2+\[Delta]r2))/(r1+\[Delta]r1+r2+\[Delta]r2);
	xTilde=Sqrt[1-(z1+\[Delta]z1)^2];
	orbitTilde = KerrGeodesics`KerrGeoOrbit`KerrGeoOrbit[a, pTilde, eTilde, xTilde, "Parametrization" -> "Phases", "Method" -> "Analytic"];
	\[Tau]Tilde = ProperTime[orbitTilde];
	frequencies = orbitTilde["Frequencies"];
	\[CapitalUpsilon]r = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"];
	\[CapitalUpsilon]z = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"];
	\[CapitalUpsilon]\[Phi] = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"];
	\[CapitalUpsilon]t = frequencies["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"] - 3*spar/(2*sqrtK)*\[Tau]Tilde["ProperTimeFrequency"];
	\[CapitalUpsilon]tz = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencyt\[Theta][a,pTilde,eTilde,xTilde,{EnTilde,JzTilde,KTilde-(JzTilde-a*EnTilde)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,pTilde,eTilde,xTilde]] + a*JzTilde - 3*spar/(2*sqrtK)*\[Tau]Tilde["\[CapitalUpsilon]\[Tau]z"];
	\[CapitalUpsilon]\[Phi]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Phi]\[Theta][a,pTilde,eTilde,xTilde,{EnTilde,JzTilde,KTilde-(JzTilde-a*EnTilde)^2},KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a,pTilde,eTilde,xTilde]] - a*EnTilde;
	\[CapitalDelta]tr = Function[qr, orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]tr"][qr] - 3*spar/(2*sqrtK)*\[Tau]Tilde["ProperTime"][0,qr,0]];
	\[CapitalDelta]tz = Function[qz, orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"][qz] - 3*spar/(2*sqrtK)*\[Tau]Tilde["ProperTime"][0,0,qz]];
	rTilde = orbitTilde["Trajectory"][[2]];
	zTilde[qz_] := Cos[orbitTilde["Trajectory"][[3]][qz]];
	UrTilde[qr_]:=If[Mod[qr,2Pi]<Pi,1,-1]*Sqrt[((rTilde[qr]^2+a^2)*EnTilde - a*JzTilde)^2 - (rTilde[qr]^2 - 2*rTilde[qr] + a^2)*(KTilde + rTilde[qr]^2)];
	UzTilde[qz_]:=If[Mod[qz,2Pi]<Pi,-1,1]*Sqrt[-((1-zTilde[qz]^2)*a*EnTilde - JzTilde)^2 + (1-zTilde[qz]^2)*(KTilde - a^2*zTilde[qz]^2)];
	Trajectoryg = {Function[{qt,qr,qz},orbitTilde["Trajectory"][[1]][qt,qr,qz] - 3*spar/(2*sqrtK)*\[Tau]Tilde["ProperTime"][0,qr,qz]],
	rTilde,zTilde,orbitTilde["Trajectory"][[4]]};
	\[Delta]x = Function[{qr,qz},{
		\[Delta]t3[rTilde[qr], zTilde[qz], UrTilde[qr], UzTilde[qz], KTilde, a],
		\[Delta]r3[rTilde[qr], zTilde[qz], EnTilde, JzTilde, KTilde, a],
		\[Delta]z3[rTilde[qr], zTilde[qz], EnTilde, JzTilde, KTilde, a],
		\[Delta]\[Phi]3[rTilde[qr], zTilde[qz], UrTilde[qr], UzTilde[qz], KTilde, a]}];
	<|
		"a"->a,
		"p"->p,
		"e"->e,
		"x"->x,
		"spar"->spar,
		"Energy"->En,
		"AngularMomentum"->Jz,
		"CarterConstantK"->K,
		"EnTilde"->EnTilde,
		"JzTilde"->JzTilde,
		"KTilde"->KTilde,
		"Frequencies"-><|
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"->\[CapitalUpsilon]r,
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(z\)]\)"->\[CapitalUpsilon]z,
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"->\[CapitalUpsilon]\[Phi],
			"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"->\[CapitalUpsilon]t
		|>,
		"ProperTimeFrequency"->\[Tau]Tilde["ProperTimeFrequency"],
		"BLFrequencies"->{\[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]}/\[CapitalUpsilon]t,
		"TrajectoryDeltas"-><|
			"\[CapitalDelta]tr"->\[CapitalDelta]tr,
			"\[CapitalDelta]tz"->\[CapitalDelta]tz,
			"\[CapitalDelta]\[Phi]r"->orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"],
			"\[CapitalDelta]\[Phi]z"->orbitTilde["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"]
		|>,
		"Trajectoryg"->Trajectoryg,
		"\[Delta]x"->\[Delta]x,
		"FourVelocitiesg"->{Function[qr,UrTilde[qr]],Function[qz,UzTilde[qz]]},
		"rootsTilde"->{r2+\[Delta]r2,r1+\[Delta]r1},
		"\[CapitalUpsilon]tz"->\[CapitalUpsilon]tz,
		"\[CapitalUpsilon]\[Phi]z"->\[CapitalUpsilon]\[Phi]z
	|>
]*)


(* ::Subsection::Closed:: *)
(*Numerical generic orbit*)


Correction[orbitGeo_,nmax_,kmax_?EvenQ]:=Module[{a,p,e,x,\[ScriptCapitalI],\[ScriptCapitalR]tfunc,\[ScriptCapitalR]\[Phi]func,\[ScriptCapitalF]rfunc,\[ScriptCapitalF]\[GothicR]func,\[ScriptCapitalG]rfunc,\[ScriptCapitalG]\[GothicR]func,\[ScriptCapitalG]\[Theta]func,\[ScriptCapitalG]\[GothicZ]func,\[ScriptCapitalH]rfunc,\[ScriptCapitalH]\[GothicR]func,\[ScriptCapitalH]\[Theta]func,\[ScriptCapitalH]\[GothicZ]func,
    \[ScriptCapitalI]1rfunc,\[ScriptCapitalI]1\[Theta]func,\[ScriptCapitalI]2func,\[ScriptCapitalI]3func,\[ScriptCapitalJ]func,\[ScriptCapitalQ]\[Theta]func,\[ScriptCapitalQ]\[GothicZ]func,\[ScriptCapitalS]rfunc,\[ScriptCapitalS]\[GothicR]func,\[ScriptCapitalS]\[Theta]func,\[ScriptCapitalS]\[GothicZ]func,\[ScriptCapitalT]rfunc,\[ScriptCapitalT]\[GothicR]func,\[ScriptCapitalT]\[Theta]func,\[ScriptCapitalT]\[GothicZ]func,\[ScriptCapitalU]1rfunc,\[ScriptCapitalU]1\[Theta]func,\[ScriptCapitalU]2func,\[ScriptCapitalU]3func,\[ScriptCapitalV]func,
    \[ScriptCapitalK]rfunc,\[ScriptCapitalK]\[GothicR]func,\[ScriptCapitalK]\[Theta]func,\[ScriptCapitalK]\[GothicZ]func,\[ScriptCapitalM]rfunc,\[ScriptCapitalM]\[GothicR]func,\[ScriptCapitalM]\[Theta]func,\[ScriptCapitalM]\[GothicZ]func,\[ScriptCapitalN]1rfunc,\[ScriptCapitalN]1\[Theta]func,\[ScriptCapitalN]2func,\[ScriptCapitalN]3func,\[ScriptCapitalP]func,Ehat,Lzhat,Qhat,Khat,\[CapitalUpsilon]rhat,\[CapitalUpsilon]\[Theta]hat,rhat,zhat,drd\[Lambda],dzd\[Lambda],
    \[CapitalDelta],\[CapitalSigma],d2rd\[Lambda]2,d2zd\[Lambda]2,R,Z,\[Chi]rhat,\[Chi]\[Theta]hat,\[Delta]udst,\[Delta]uds\[Phi],\[CapitalLambda]rhat,\[CapitalLambda]\[Theta]hat,stepsr,steps\[Theta],ExpniTable,ExpjkTable,\[ScriptCapitalR]tlist,\[ScriptCapitalR]tfour,\[ScriptCapitalR]\[Phi]list,\[ScriptCapitalR]\[Phi]four,\[ScriptCapitalF]rfour,\[ScriptCapitalF]rlist,\[ScriptCapitalF]\[GothicR]four,\[ScriptCapitalF]\[GothicR]list,
    \[ScriptCapitalG]rfour,\[ScriptCapitalG]rlist,\[ScriptCapitalG]\[GothicR]four,\[ScriptCapitalG]\[GothicR]list,\[ScriptCapitalG]\[Theta]four,\[ScriptCapitalG]\[Theta]list,\[ScriptCapitalG]\[GothicZ]four,\[ScriptCapitalG]\[GothicZ]list,\[ScriptCapitalH]rfour,\[ScriptCapitalH]rlist,\[ScriptCapitalH]\[GothicR]four,\[ScriptCapitalH]\[GothicR]list,\[ScriptCapitalH]\[Theta]four,\[ScriptCapitalH]\[Theta]list,\[ScriptCapitalH]\[GothicZ]four,\[ScriptCapitalH]\[GothicZ]list,\[ScriptCapitalI]1rfour,\[ScriptCapitalI]1rlist,\[ScriptCapitalI]1\[Theta]four,\[ScriptCapitalI]1\[Theta]list,\[ScriptCapitalI]2four,\[ScriptCapitalI]2list,
    \[ScriptCapitalI]3four,\[ScriptCapitalI]3list,\[ScriptCapitalJ]four,\[ScriptCapitalJ]list,\[ScriptCapitalQ]\[Theta]four,\[ScriptCapitalQ]\[Theta]list,\[ScriptCapitalQ]\[GothicZ]four,\[ScriptCapitalQ]\[GothicZ]list,\[ScriptCapitalS]rfour,\[ScriptCapitalS]rlist,\[ScriptCapitalS]\[GothicR]four,\[ScriptCapitalS]\[GothicR]list,\[ScriptCapitalS]\[Theta]four,\[ScriptCapitalS]\[Theta]list,\[ScriptCapitalS]\[GothicZ]four,\[ScriptCapitalS]\[GothicZ]list,\[ScriptCapitalT]rfour,\[ScriptCapitalT]rlist,\[ScriptCapitalT]\[GothicR]four,\[ScriptCapitalT]\[GothicR]list,\[ScriptCapitalT]\[Theta]four,\[ScriptCapitalT]\[Theta]list,
    \[ScriptCapitalT]\[GothicZ]four,\[ScriptCapitalT]\[GothicZ]list,\[ScriptCapitalU]1rfour,\[ScriptCapitalU]1rlist,\[ScriptCapitalU]1\[Theta]four,\[ScriptCapitalU]1\[Theta]list,\[ScriptCapitalU]2four,\[ScriptCapitalU]2list,\[ScriptCapitalU]3four,\[ScriptCapitalU]3list,\[ScriptCapitalV]four,\[ScriptCapitalV]list,\[ScriptCapitalK]rfour,\[ScriptCapitalK]rlist,\[ScriptCapitalK]\[GothicR]four,\[ScriptCapitalK]\[GothicR]list,\[ScriptCapitalK]\[Theta]four,\[ScriptCapitalK]\[Theta]list,\[ScriptCapitalK]\[GothicZ]four,\[ScriptCapitalK]\[GothicZ]list,\[ScriptCapitalM]rfour,\[ScriptCapitalM]rlist,
    \[ScriptCapitalM]\[GothicR]four,\[ScriptCapitalM]\[GothicR]list,\[ScriptCapitalM]\[Theta]four,\[ScriptCapitalM]\[Theta]list,\[ScriptCapitalM]\[GothicZ]four,\[ScriptCapitalM]\[GothicZ]list,\[ScriptCapitalN]1rfour,\[ScriptCapitalN]1rlist,\[ScriptCapitalN]1\[Theta]four,\[ScriptCapitalN]1\[Theta]list,\[ScriptCapitalN]2four,\[ScriptCapitalN]2list,\[ScriptCapitalN]3four,\[ScriptCapitalN]3list,\[ScriptCapitalP]four,\[ScriptCapitalP]list,M1r,M1\[Theta],M1\[GothicR],M1\[GothicZ],M1\[CapitalUpsilon]r,M1\[CapitalUpsilon]\[Theta],M1t,M1\[Phi],M1,
    M2r,M2\[Theta],M2\[GothicR],M2\[GothicZ],M2\[CapitalUpsilon]r,M2\[CapitalUpsilon]\[Theta],M2t,M2\[Phi],M2,M3r,M3\[Theta],M3\[GothicR],M3\[GothicZ],M3\[CapitalUpsilon]r,M3\[CapitalUpsilon]\[Theta],M3t,M3\[Phi],M3,M,c,v,\[Delta]\[Chi]rSfourIm,\[Delta]\[Chi]\[Theta]SfourIm,\[Delta]\[GothicR]four,\[Delta]\[GothicZ]four,\[CapitalUpsilon]rS,\[CapitalUpsilon]\[Theta]S,udt0S,ud\[Phi]0S},
  drd\[Lambda][\[Lambda]r_]:=Sign[Pi/\[CapitalUpsilon]rhat-Mod[\[Lambda]r,2Pi/\[CapitalUpsilon]rhat]]*Sqrt[(Ehat*(rhat[\[Lambda]r]^2+a^2)-a*Lzhat)^2-\[CapitalDelta][rhat[\[Lambda]r]]*(rhat[\[Lambda]r]^2+Khat)];(* geodesic radial velocity *)
  dzd\[Lambda][\[Lambda]\[Theta]_]:=-Sign[Pi/\[CapitalUpsilon]\[Theta]hat-Mod[\[Lambda]\[Theta],2Pi/\[CapitalUpsilon]\[Theta]hat]]*Sqrt[Qhat-zhat[\[Lambda]\[Theta]]^2*(a^2*(1-Ehat^2)*(1-zhat[\[Lambda]\[Theta]]^2)+Lzhat^2+Qhat)];(* geodesic polar velocity *)
  \[CapitalDelta][r_]:=r^2-2*r+a^2;
  \[CapitalSigma][r_,z_]:=r^2+a^2*z^2;
  d2rd\[Lambda]2[r_]:=1/2 (-((-2+2 r) (Khat+r^2))-2 r (a^2-2 r+r^2)+4 Ehat r (-a Lzhat+Ehat (a^2+r^2)));(* second derivative of r wrt \[Lambda] *)
  d2zd\[Lambda]2[z_]:=(2 a^2 (1-Ehat^2) z^3-2 z (Lzhat^2+Qhat+a^2 (1-Ehat^2) (1-z^2)))/2;(* second derivative of z=cos(\[Theta]) wrt \[Lambda] *)
  R[r_]:=((r^2+a^2)*Ehat-a*Lzhat)^2-\[CapitalDelta][r]*(Khat+r^2);(* equations of motion for r,z *)
  Z[z_]:=Qhat-z^2*(a^2*(1-Ehat^2)*(1-z^2)+Lzhat^2+Qhat);
  \[Chi]rhat[\[Lambda]r_]:=\[Chi]rhat[\[Lambda]r]=If[Mod[\[Lambda]r,2Pi/\[CapitalUpsilon]rhat]<Pi/\[CapitalUpsilon]rhat,ArcCos[(p/orbitGeo["Trajectory"][[2]][\[Lambda]r]-1)/e]+2Pi*Floor[\[Lambda]r/(2Pi/\[CapitalUpsilon]rhat)],2Pi+2Pi*Floor[\[Lambda]r/(2Pi/\[CapitalUpsilon]rhat)]-ArcCos[(p/orbitGeo["Trajectory"][[2]][\[Lambda]r]-1)/e]];(* calculation of \[Chi]rhat from rhat *)
  \[Chi]\[Theta]hat[\[Lambda]\[Theta]_]:=\[Chi]\[Theta]hat[\[Lambda]\[Theta]]=If[Mod[\[Lambda]\[Theta],2Pi/\[CapitalUpsilon]\[Theta]hat]<Pi/\[CapitalUpsilon]\[Theta]hat,ArcCos[Cos[orbitGeo["Trajectory"][[3]][\[Lambda]\[Theta]]]/Sin[\[ScriptCapitalI]]]+2Pi*Floor[\[Lambda]\[Theta]/(2Pi/\[CapitalUpsilon]\[Theta]hat)],2Pi+2Pi*Floor[\[Lambda]\[Theta]/(2Pi/\[CapitalUpsilon]\[Theta]hat)]-ArcCos[Cos[orbitGeo["Trajectory"][[3]][\[Lambda]\[Theta]]]/Sin[\[ScriptCapitalI]]]];
  \[Delta]udst[\[Lambda]r_,\[Lambda]\[Theta]_]:=\[Delta]udst[\[Lambda]r,\[Lambda]\[Theta]]=Re[Sum[If[k==0&&n==0,0,Im[\[ScriptCapitalR]tfour[[n+nmax+1,k+kmax+1]]]/(-k*\[CapitalUpsilon]\[Theta]hat-n*\[CapitalUpsilon]rhat)*Exp[-I*k*\[CapitalUpsilon]\[Theta]hat*\[Lambda]\[Theta]-I*n*\[CapitalUpsilon]rhat*\[Lambda]r]],{n,-nmax,nmax},{k,-kmax,kmax,2}]];(* definition of \[Delta]udst summed over even k for aligned spin *)
  \[Delta]uds\[Phi][\[Lambda]r_,\[Lambda]\[Theta]_]:=\[Delta]uds\[Phi][\[Lambda]r,\[Lambda]\[Theta]]=Re[Sum[If[k==0&&n==0,0,Im[\[ScriptCapitalR]\[Phi]four[[n+nmax+1,k+kmax+1]]]/(-k*\[CapitalUpsilon]\[Theta]hat-n*\[CapitalUpsilon]rhat)*Exp[-I*k*\[CapitalUpsilon]\[Theta]hat*\[Lambda]\[Theta]-I*n*\[CapitalUpsilon]rhat*\[Lambda]r]],{n,-nmax,nmax},{k,-kmax,kmax,2}]];
  \[ScriptCapitalR]tfunc[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=-((2 a^2 dzd\[Lambda] r z (3 r^2 (2 Khat+r^2)-2 a^2 (3 Khat+4 r^2) z^2+a^4 z^4)+drd\[Lambda] (3 Khat r^4-9 a^2 r^2 (2 Khat+r^2) z^2+a^4 (3 Khat+14 r^2) z^4-a^6 z^6))/(Sqrt[Khat] (r^2+a^2 z^2)^4));(* functionf for t and \[Phi] equations *)
  \[ScriptCapitalR]\[Phi]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=-((a (-2 dzd\[Lambda] r (a^2+r^2) z (3 r^2 (2 Khat+r^2)-2 a^2 (3 Khat+4 r^2) z^2+a^4 z^4)+drd\[Lambda] (-1+z^2) (3 Khat r^4-9 a^2 r^2 (2 Khat+r^2) z^2+a^4 (3 Khat+14 r^2) z^4-a^6 z^6)))/(Sqrt[Khat] (r^2+a^2 z^2)^4));
  \[ScriptCapitalF]rfunc[\[Chi]rhat_]:=(e p Sin[\[Chi]rhat])/(1+e Cos[\[Chi]rhat])^2;(* functions for r equations *)
  \[ScriptCapitalF]\[GothicR]func[]:=1;
  \[ScriptCapitalG]rfunc[r_,z_,drd\[Lambda]_,\[Chi]rhat_]:=(2 e r (r^3 (p (-a^2+p)-(a^2 (-1+e^2)+p^2) r+(-1+e^2+p) r^2)+a^2 (a^2 p^2+r (-p (3 a^2+p)+r (-2 a^2 (-1+e^2)+4 p-r (3-3 e^2+p+(-1+e^2) r)))) z^2) drd\[Lambda] Sin[\[Chi]rhat])/(p (a^2+(-2+r) r) (p^2+r (-2 p+r-e^2 r)) (r^2+a^2 z^2));
  \[ScriptCapitalG]\[GothicR]func[r_,z_,drd\[Lambda]_]:=(-((2 (-1+r))/(a^2+(-2+r) r))-(2 r)/(r^2+a^2 z^2))*drd\[Lambda];
  \[ScriptCapitalG]\[Theta]func[r_,z_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=-((2 dzd\[Lambda] r (a^2+(-2+r) r) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/((-1+z^2) (r^2+a^2 z^2)));
  \[ScriptCapitalG]\[GothicZ]func[r_,z_,dzd\[Lambda]_]:=(2 r (a^2+(-2+r) r) dzd\[Lambda])/((-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalH]rfunc[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,d2rd\[Lambda]2_,\[Chi]rhat_]:=1/(p (a^2+(-2+r) r)^2 (-1+z^2) (r^2+a^2 z^2)^2) e r^2 (Lzhat^2 (-2+r)^2 r^6+4 a Ehat Lzhat r^6 (-1+z^2)+2 Ehat^2 (-3+r) r^8 (-1+z^2)-4 a^5 Ehat Lzhat r^2 z^2 (-10+8 r-z^2) (-1+z^2)-4 a^3 Ehat Lzhat r^3 (-1+z^2) (r (-1+2 r)+2 (-4+r) (-1+r) z^2)+a^2 r^3 (4 Ehat^2 r^3 (-1+z^2) (r+(-3+r) z^2)+Lzhat^2 (r (2+3 (-2+r) r)+(-16+r (22+(-10+r) r)) z^2))+4 a^7 Ehat Lzhat z^2 (-z^2+z^4-2 r (-1+z^4))+a^6 z^2 (-4 Ehat^2 r^2 (-1+z^2) (5-3 r+(-4+r) z^2)+Lzhat^2 (-4 r+(2+r (-6+5 r)) z^2-(2+(-2+r) r) z^4))+a^8 z^2 (Lzhat^2 z^4-2 Ehat^2 (-1+z^2) (z^2+r (-2+z^2)))+a^4 r^2 (-2 Ehat^2 r (-1+z^2) (-((-1+r) r)-2 (4+r (-4+3 r)) z^2+(8+(-5+r) r) z^4)+
    Lzhat^2 z^2 (20-r (20+r (-7+z^2))))) Sin[\[Chi]rhat]+1/(e p r) Csc[\[Chi]rhat] ((-p^2+3 p r+2 (-1+e^2) r^2) d2rd\[Lambda]2+(-2 (-1+e^2+p)-(-1+e^2) r-(e p^2)/(2 (p+(-1+e) r))+(e p^2)/(2 p-2 (1+e) r)+(2 (-1+a^2) (2 a^2 (-1+e^2+p)+((-4+a^2) (-1+e^2)-4 p+p^2) r))/(a^2+(-2+r) r)^2+(((-2+a^2) (-1+e^2)-2 p+p^2) (-2+r))/(a^2+(-2+r) r)+(2 r^3 (p+(-1+e) r) (-p+(1+e) r))/(r^2+a^2 z^2)^2+(r (3 p^2-8 p r-5 (-1+e^2) r^2))/(r^2+a^2 z^2)) drd\[Lambda]^2)+(e r^2 (-a^2 r^2+r^4+a^2 (a^2+r (-4+3 r)) z^2) dzd\[Lambda]^2 Sin[\[Chi]rhat])/(p (-1+z^2) (r^2+a^2 z^2)^2);
  \[ScriptCapitalH]\[GothicR]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=1/2 (1/(a^2+(-2+r) r)^2 (2 a^2 (-2 a^2 Ehat^2+4 a Ehat Lzhat+(-2+a^2) Lzhat^2)+2 r (2 a^2 (-a^2 Ehat^2-4 a Ehat Lzhat+Lzhat^2)+r (a (16 a Ehat^2+4 Ehat Lzhat-a Lzhat^2)-2 Ehat^2 r (8+2 a^2-5 r+r^2))))+(Lzhat^2 (a^4-4 a^2 r+2 a^2 r^2+r^4))/((a^2+r^2)^2 (-1+z))-(Lzhat^2 (a^4-4 a^2 r+2 a^2 r^2+r^4))/((a^2+r^2)^2 (1+z))-(8 r^3 (a (a Ehat-Lzhat)+Ehat r^2)^2)/((a^2+r^2) (r^2+a^2 z^2)^2)+(8 r (a (a Ehat-Lzhat)+Ehat r^2) (a^3 (a Ehat-Lzhat)+3 a^2 Ehat r^2+2 Ehat r^4))/((a^2+r^2)^2 (r^2+a^2 z^2)))+drd\[Lambda]^2 ((2-a^2+(-2+r) r)/(a^2+(-2+r) r)^2+(2 r^2)/(r^2+a^2 z^2)^2-1/(r^2+a^2 z^2))+(dzd\[Lambda]^2 (-a^2 r^2+r^4+a^2 (a^2+r (-4+3 r)) z^2))/((-1+z^2) (r^2+a^2 z^2)^2);
  \[ScriptCapitalH]\[Theta]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=1/((-1+z^2)^2 (r^2+a^2 z^2)^2) 2 r z (Lzhat^2 r^4+2 a^2 Lzhat^2 r^2 z^2+a^4 Lzhat^2 z^4+2 r^3 (-a Ehat-Lzhat+a Ehat z^2) (-a Ehat+Lzhat+a Ehat z^2)+2 a^2 r (-a Ehat+Lzhat+a Ehat z^2) (-a Ehat+Lzhat+(a Ehat-2 Lzhat) z^2)) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat]-(2 a^2 drd\[Lambda]^2 r z Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/(r^2+a^2 z^2)^2-(2 dzd\[Lambda]^2 r (a^2+(-2+r) r) z (r^2 Cot[\[ScriptCapitalI]]^2+a^2 (1-2 z^2+z^4 Csc[\[ScriptCapitalI]]^2)) Csc[\[Chi]\[Theta]hat] Sin[\[ScriptCapitalI]])/((-1+z^2)^2 (r^2+a^2 z^2)^2);
  \[ScriptCapitalH]\[GothicZ]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=(2 a^2 drd\[Lambda]^2 r z)/(r^2+a^2 z^2)^2-(2 dzd\[Lambda]^2 r (a^2+(-2+r) r) z (r^2+a^2 (-1+2 z^2)))/((-1+z^2)^2 (r^2+a^2 z^2)^2)+1/((-1+z^2)^2 (r^2+a^2 z^2)^2) 2 r z (-Lzhat^2 r^4-2 a^2 Lzhat^2 r^2 z^2-a^4 Lzhat^2 z^4-2 r^3 (-a Ehat-Lzhat+a Ehat z^2) (-a Ehat+Lzhat+a Ehat z^2)-2 a^2 r (-a Ehat+Lzhat+a Ehat z^2) (-a Ehat+Lzhat+(a Ehat-2 Lzhat) z^2));
  \[ScriptCapitalI]1rfunc[r_,z_,drd\[Lambda]_,d2rd\[Lambda]2_]:=(2 (r^2 (drd\[Lambda]^2 (3-2 r)+d2rd\[Lambda]2 (-2+r) r)+a^4 d2rd\[Lambda]2 z^2+a^2 (d2rd\[Lambda]2 r (r+(-2+r) z^2)-drd\[Lambda]^2 (r+(-1+r) z^2))))/((a^2+(-2+r) r) (r^2+a^2 z^2) \[CapitalUpsilon]rhat);
  \[ScriptCapitalI]1\[Theta]func[r_,z_,dzd\[Lambda]_]:=(2 dzd\[Lambda]^2 r (a^2+(-2+r) r)/\[CapitalUpsilon]\[Theta]hat)/((-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalI]2func[r_,z_]:=1/(a^2+(-2+r) r) 2 (a^3 (a Ehat-Lzhat)+r^2 (a (2 a Ehat+Lzhat)+Ehat (-4+r) r))-(4 r^2 (a (a Ehat-Lzhat)+Ehat r^2))/(r^2+a^2 z^2);
  \[ScriptCapitalI]3func[r_,z_]:=(2 (Lzhat (-2+r)^2 r^3-a Ehat r^3 (-4+3 r) (-1+z^2)+a^5 Ehat z^2 (-1+z^2)+a^4 Lzhat z^2 (1+(-1+r) z^2)+a^2 Lzhat r^2 (-1+(-3+2 r) z^2)-a^3 Ehat r^2 (-1+z^4)))/((a^2+(-2+r) r) (-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalJ]func[r_,z_,\[Delta]udst_,\[Delta]uds\[Phi]_]:=((a^2 Ehat-a Lzhat+Ehat r^2) (-3 Khat r^4+9 a^2 r^2 (2 Khat+r^2) z^2-a^4 (3 Khat+14 r^2) z^4+a^6 z^6))/(Sqrt[Khat] (r^2+a^2 z^2)^3)+\[ScriptCapitalI]2func[r,z]*\[Delta]udst+\[ScriptCapitalI]3func[r,z]*\[Delta]uds\[Phi];
  \[ScriptCapitalK]rfunc[r_,z_,drd\[Lambda]_,\[Chi]rhat_]:=(2 drd\[Lambda] e r^2 Sin[\[Chi]rhat])/(p (a^2+(-2+r) r) (r^2+a^2 z^2));(* functions for norm equation *)
  \[ScriptCapitalK]\[GothicR]func[r_,z_,drd\[Lambda]_]:=(2 drd\[Lambda])/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalK]\[Theta]func[r_,z_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=(2 dzd\[Lambda] Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/((-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalK]\[GothicZ]func[r_,z_,dzd\[Lambda]_]:=-((2 dzd\[Lambda])/((-1+z^2) (r^2+a^2 z^2)));
  \[ScriptCapitalM]rfunc[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,\[Chi]rhat_]:=(2 drd\[Lambda]^2 (r^3 ((a^2-p) p+r (a^2 (-1+e^2)+p^2-(-1+e^2+p) r))+a^2 (-a^2 p^2+r (p (3 a^2+p)+r (2 a^2 (-1+e^2)-4 p+r (3-3 e^2+p+(-1+e^2) r)))) z^2) Csc[\[Chi]rhat])/(e p r (a^2+(-2+r) r)^2 (r^2+a^2 z^2)^2)+(2 dzd\[Lambda]^2 e r^3 Sin[\[Chi]rhat])/(p (-1+z^2) (r^2+a^2 z^2)^2)+1/(p (a^2+(-2+r) r)^2 (-1+z^2) (r^2+a^2 z^2)^2) 2 e r^2 (Lzhat^2 r^5+a^4 Lzhat^2 r z^4+Ehat^2 r^6 (-1+z^2)-a^4 (-a Ehat+Lzhat)^2 z^2 (-1+z^2)-r^4 (2 (a Ehat-2 Lzhat) (a Ehat-Lzhat)+a Ehat z^2 (-3 a Ehat+6 Lzhat+a Ehat z^2))+a^2 r^2 (-(-a Ehat+Lzhat)^2-(a Ehat+Lzhat) z^2 (-3 a Ehat+3 Lzhat+2 a Ehat z^2))+2 r^3 (2 (-a Ehat+Lzhat)^2+a z^2 (-4 a Ehat^2+4 Ehat Lzhat+a Lzhat^2+2 a Ehat^2 z^2))) Sin[\[Chi]rhat];
  \[ScriptCapitalM]\[GothicR]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=(2 dzd\[Lambda]^2 r)/((-1+z^2) (r^2+a^2 z^2)^2)+(2 drd\[Lambda]^2 (-r (a^2+r (-3+2 r))-a^2 (-1+r) z^2))/((a^2+(-2+r) r)^2 (r^2+a^2 z^2)^2)+1/((a^2+(-2+r) r)^2 (-1+z^2) (r^2+a^2 z^2)^2) 2 (Lzhat^2 r^5+a^4 Lzhat^2 r z^4+Ehat^2 r^6 (-1+z^2)-a^4 (-a Ehat+Lzhat)^2 z^2 (-1+z^2)-r^4 (2 (a Ehat-2 Lzhat) (a Ehat-Lzhat)+a Ehat z^2 (-3 a Ehat+6 Lzhat+a Ehat z^2))+a^2 r^2 (-(-a Ehat+Lzhat)^2-(a Ehat+Lzhat) z^2 (-3 a Ehat+3 Lzhat+2 a Ehat z^2))+2 r^3 (2 (-a Ehat+Lzhat)^2+a z^2 (-4 a Ehat^2+4 Ehat Lzhat+a Lzhat^2+2 a Ehat^2 z^2)));
  \[ScriptCapitalM]\[Theta]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=(2 a^2 drd\[Lambda]^2 z Sin[\[Chi]\[Theta]hat] Sin[\[ScriptCapitalI]])/((a^2+(-2+r) r) (r^2+a^2 z^2)^2)-1/((a^2+(-2+r) r) (-1+z^2)^2 (r^2+a^2 z^2)^2) 2 z (Lzhat^2 r^4+2 a^2 Lzhat^2 r^2 z^2+a^4 Lzhat^2 z^4+2 r^3 (-a Ehat-Lzhat+a Ehat z^2) (-a Ehat+Lzhat+a Ehat z^2)+2 a^2 r (-a Ehat+Lzhat+a Ehat z^2) (-a Ehat+Lzhat+(a Ehat-2 Lzhat) z^2)) Sin[\[Chi]\[Theta]hat] Sin[\[ScriptCapitalI]]+(2 dzd\[Lambda]^2 z (r^2 Cot[\[ScriptCapitalI]]^2+a^2 (1-2 z^2+z^4 Csc[\[ScriptCapitalI]]^2)) Sin[\[ScriptCapitalI]])/((-1+z^2)^2 (r^2+a^2 z^2)^2 Sin[\[Chi]\[Theta]hat]);
  \[ScriptCapitalM]\[GothicZ]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=-((2 a^2 drd\[Lambda]^2 z)/((a^2+(-2+r) r) (r^2+a^2 z^2)^2))+(2 dzd\[Lambda]^2 z (r^2+a^2 (-1+2 z^2)))/((-1+z^2)^2 (r^2+a^2 z^2)^2)+1/((a^2+(-2+r) r) (-1+z^2)^2 (r^2+a^2 z^2)^2) 2 z (Lzhat^2 r^4+2 a^2 Lzhat^2 r^2 z^2+a^4 Lzhat^2 z^4+2 r^3 (-a Ehat-Lzhat+a Ehat z^2) (-a Ehat+Lzhat+a Ehat z^2)+2 a^2 r (-a Ehat+Lzhat+a Ehat z^2) (-a Ehat+Lzhat+(a Ehat-2 Lzhat) z^2));
  \[ScriptCapitalN]1rfunc[r_,z_,drd\[Lambda]_]:=(2 drd\[Lambda]^2/\[CapitalUpsilon]rhat)/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalN]1\[Theta]func[r_,z_,dzd\[Lambda]_]:=-((2 dzd\[Lambda]^2/\[CapitalUpsilon]\[Theta]hat)/((-1+z^2) (r^2+a^2 z^2)));
  \[ScriptCapitalN]2func[r_,z_]:=2 Ehat+(4 r (a (a Ehat-Lzhat)+Ehat r^2))/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalN]3func[r_,z_]:=(-2 r (2 a Ehat-2 Lzhat+Lzhat r)-2 a (a Lzhat-2 Ehat r) z^2)/((a^2+(-2+r) r) (-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalP]func[r_,z_,\[Delta]udst_,\[Delta]uds\[Phi]_]:=\[ScriptCapitalN]2func[r,z]*\[Delta]udst+\[ScriptCapitalN]3func[r,z]*\[Delta]uds\[Phi];
  \[ScriptCapitalQ]\[Theta]func[z_,\[Chi]\[Theta]hat_]:=(Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/Sqrt[1-z^2];(* functions for \[Theta] equation *)
  \[ScriptCapitalQ]\[GothicZ]func[z_]:=-(1/Sqrt[1-z^2]);
  \[ScriptCapitalS]rfunc[r_,z_,drd\[Lambda]_,\[Chi]rhat_]:=(2 a^2 e r^2 z Sqrt[1-z^2] drd\[Lambda] Sin[\[Chi]rhat])/(p (a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalS]\[GothicR]func[r_,z_,drd\[Lambda]_]:=(2 a^2 z Sqrt[1-z^2] drd\[Lambda])/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalS]\[Theta]func[r_,z_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=-((2 z Sqrt[1-z^2] (r^2 Cot[\[ScriptCapitalI]]^2+a^2 (1-2 z^2+z^4 Csc[\[ScriptCapitalI]]^2)) dzd\[Lambda] Sin[\[ScriptCapitalI]])/((-1+z^2)^2 (r^2+a^2 z^2) Sin[\[Chi]\[Theta]hat]));
  \[ScriptCapitalS]\[GothicZ]func[r_,z_,dzd\[Lambda]_]:=-((2 dzd\[Lambda] z (r^2+a^2 (-1+2 z^2)))/((1-z^2)^(3/2) (r^2+a^2 z^2)));
  \[ScriptCapitalT]rfunc[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,\[Chi]rhat_]:=1/(p (a^2+(-2+r) r)^2 Sqrt[1-z^2] (r^2+a^2 z^2)^2) 2 a^2 e r^2 z (-Lzhat^2 r^5-a^4 Lzhat^2 r z^4-Ehat^2 r^6 (-1+z^2)+a^4 (-a Ehat+Lzhat)^2 z^2 (-1+z^2)+r^4 (2 (a Ehat-2 Lzhat) (a Ehat-Lzhat)+a Ehat z^2 (-3 a Ehat+6 Lzhat+a Ehat z^2))+a^2 r^2 ((-a Ehat+Lzhat)^2+(a Ehat+Lzhat) z^2 (-3 a Ehat+3 Lzhat+2 a Ehat z^2))+2 r^3 (-2 (-a Ehat+Lzhat)^2-a z^2 (-4 a Ehat^2+4 Ehat Lzhat+a Lzhat^2+2 a Ehat^2 z^2))) Sin[\[Chi]rhat]+(2 a^2 drd\[Lambda]^2 e r z Sqrt[1-z^2] (r^3 (p (-a^2+p)-(a^2 (-1+e^2)+p^2) r+(-1+e^2+p) r^2)+a^2 (a^2 p^2+r (-p (3 a^2+p)+r (-2 a^2 (-1+e^2)+4 p-r (3-3 e^2+p+(-1+e^2) r)))) z^2) Sin[\[Chi]rhat])/((a^2+(-2+r) r)^2 (p^3+p r (-2 p+r-e^2 r)) (r^2+a^2 z^2)^2)
    -(2 a^2 dzd\[Lambda]^2 e r^3 z Sin[\[Chi]rhat])/(p Sqrt[1-z^2] (r^2+a^2 z^2)^2);
  \[ScriptCapitalT]\[GothicR]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=1/((a^2+(-2+r) r)^2 Sqrt[1-z^2] (r^2+a^2 z^2)^2) 2 a^2 z (-Lzhat^2 r^5-a^4 Lzhat^2 r z^4-Ehat^2 r^6 (-1+z^2)+a^4 (-a Ehat+Lzhat)^2 z^2 (-1+z^2)+r^4 (2 (a Ehat-2 Lzhat) (a Ehat-Lzhat)+a Ehat z^2 (-3 a Ehat+6 Lzhat+a Ehat z^2))+a^2 r^2 ((-a Ehat+Lzhat)^2+(a Ehat+Lzhat) z^2 (-3 a Ehat+3 Lzhat+2 a Ehat z^2))+2 r^3 (-2 (-a Ehat+Lzhat)^2-a z^2 (-4 a Ehat^2+4 Ehat Lzhat+a Lzhat^2+2 a Ehat^2 z^2)))-(2 a^2 drd\[Lambda]^2 z Sqrt[1-z^2] (r (a^2+r (-3+2 r))+a^2 (-1+r) z^2))/((a^2+(-2+r) r)^2 (r^2+a^2 z^2)^2)-(2 a^2 dzd\[Lambda]^2 r z)/(Sqrt[1-z^2] (r^2+a^2 z^2)^2);
  \[ScriptCapitalT]\[Theta]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,d2zd\[Lambda]2_,\[Chi]\[Theta]hat_]:=1/((a^2+(-2+r) r) (1-z^2)^(5/2) (r^2+a^2 z^2)^2) (3 a^6 Lzhat^2 z^6+Lzhat^2 r^6 (1+2 z^2)+a^4 Lzhat^2 r^2 z^4 (7+2 z^2)+a^2 Lzhat^2 r^4 z^2 (5+4 z^2)+2 a^4 r z^2 (-(-a Ehat+Lzhat)^2+2 (-a Ehat+Lzhat)^2 z^2-(a^2 Ehat^2-2 a Ehat Lzhat+4 Lzhat^2) z^4)+2 a^2 r^3 ((-a Ehat+Lzhat)^2-(5 a Ehat-3 Lzhat) (a Ehat-Lzhat) z^2+(7 a^2 Ehat^2-10 a Ehat Lzhat-4 Lzhat^2) z^4+a Ehat (-3 a Ehat+4 Lzhat) z^6)-2 r^5 (-a^2 Ehat^2+Lzhat^2+2 (2 a^2 Ehat^2+Lzhat^2) z^2+a^2 Ehat^2 z^4 (-5+2 z^2))) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat]+(a^2 drd\[Lambda]^2 (a^2 z^2+r^2 (-1+2 z^2)) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/((a^2+(-2+r) r) Sqrt[1-z^2] (r^2+a^2 z^2)^2)
    -1/(2 (1-z^2)^(5/2) (r^2+a^2 z^2)^2) dzd\[Lambda]^2 Csc[\[ScriptCapitalI]] Csc[\[Chi]\[Theta]hat] (2 r^4 (1+2 z^2) Cos[\[ScriptCapitalI]]^2+a^2 r^2 (-4 z^6+2 z^4 (7+2 Cos[2 \[ScriptCapitalI]])+z^2 (-5+3 Cos[2 \[ScriptCapitalI]])+2 Sin[\[ScriptCapitalI]]^2)+2 a^4 z^2 (z^4 (1+2 Cos[2 \[ScriptCapitalI]])-Sin[\[ScriptCapitalI]]^2+2 z^2 Sin[\[ScriptCapitalI]]^2))+(2 z Cos[\[ScriptCapitalI]] Cot[\[ScriptCapitalI]] Csc[\[Chi]\[Theta]hat] (dzd\[Lambda]^2 z+d2zd\[Lambda]2 (-z^2+Sin[\[ScriptCapitalI]]^2)))/((1-z^2)^(3/2) (-1+2 z^2+Cos[2 \[ScriptCapitalI]]));
  \[ScriptCapitalT]\[GothicZ]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,d2zd\[Lambda]2_]:=(a^2 drd\[Lambda]^2 (r^2-(a^2+2 r^2) z^2))/((a^2+(-2+r) r) Sqrt[1-z^2] (r^2+a^2 z^2)^2)+1/((a^2+(-2+r) r) (1-z^2)^(5/2) (r^2+a^2 z^2)^2) (-3 a^6 Lzhat^2 z^6-Lzhat^2 r^6 (1+2 z^2)-a^4 Lzhat^2 r^2 z^4 (7+2 z^2)-a^2 Lzhat^2 r^4 z^2 (5+4 z^2)+2 a^4 r z^2 ((-a Ehat+Lzhat)^2-2 (-a Ehat+Lzhat)^2 z^2+(a^2 Ehat^2-2 a Ehat Lzhat+4 Lzhat^2) z^4)+2 a^2 r^3 (-(-a Ehat+Lzhat)^2+(5 a Ehat-3 Lzhat) (a Ehat-Lzhat) z^2+(-7 a^2 Ehat^2+10 a Ehat Lzhat+4 Lzhat^2) z^4+a Ehat (3 a Ehat-4 Lzhat) z^6)+2 r^5 (-a^2 Ehat^2+Lzhat^2+2 (2 a^2 Ehat^2+Lzhat^2) z^2+a^2 Ehat^2 z^4 (-5+2 z^2)))+(z (-d2zd\[Lambda]2+z (-dzd\[Lambda]^2+d2zd\[Lambda]2 z) Csc[\[ScriptCapitalI]]^2))/((1-z^2)^(3/2) (1-z^2 Csc[\[ScriptCapitalI]]^2))+(dzd\[Lambda]^2 (2 a^4 z^2 (
    -3 z^6+z^4 (3-2 Cos[2 \[ScriptCapitalI]])+z^2 (-2+Cos[2 \[ScriptCapitalI]])+Sin[\[ScriptCapitalI]]^2)+2 r^4 (-z^2 (1+z^2+Cos[2 \[ScriptCapitalI]])+Sin[\[ScriptCapitalI]]^2)-a^2 r^2 (4 z^6+z^2 (-5+3 Cos[2 \[ScriptCapitalI]])+z^4 (6+4 Cos[2 \[ScriptCapitalI]])+2 Sin[\[ScriptCapitalI]]^2)))/((1-z^2)^(5/2) (r^2+a^2 z^2)^2 (-1+2 z^2+Cos[2 \[ScriptCapitalI]]));
  \[ScriptCapitalU]1rfunc[r_,z_,drd\[Lambda]_]:=(2 a^2 z Sqrt[1-z^2] drd\[Lambda]^2/\[CapitalUpsilon]rhat)/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalU]1\[Theta]func[r_,z_,dzd\[Lambda]_,d2zd\[Lambda]2_,\[Chi]\[Theta]hat_]:=-((-2 d2zd\[Lambda]2 (-1+z^2) (r^2+a^2 z^2)+2 dzd\[Lambda]^2 z (r^2+a^2 (-1+2 z^2)))/((1-z^2)^(3/2) (r^2+a^2 z^2) \[CapitalUpsilon]\[Theta]hat));
  \[ScriptCapitalU]2func[r_,z_]:=(4 a^2 r (a^2 Ehat-a Lzhat+Ehat r^2) z Sqrt[1-z^2])/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalU]3func[r_,z_]:=-((2 z (Lzhat (-2+r) r^3+a^4 Lzhat z^4-2 a^3 Ehat r (-1+z^2)^2+2 a^2 Lzhat r (1+(-2+r) z^2)))/((a^2+(-2+r) r) (1-z^2)^(3/2) (r^2+a^2 z^2)));
  \[ScriptCapitalV]func[r_,z_,\[Delta]udst_,\[Delta]uds\[Phi]_]:=-((2 a r z (3 r^2 (2 Khat+r^2)-2 a^2 (3 Khat+4 r^2) z^2+a^4 z^4) (Lzhat+a Ehat (-1+z^2)))/(Sqrt[Khat] (r^2+a^2 z^2)^3))/Sqrt[1-z^2]+\[ScriptCapitalU]2func[r,z]*\[Delta]udst+\[ScriptCapitalU]3func[r,z]*\[Delta]uds\[Phi];
  a=orbitGeo["a"];(*orbital parameters*)
  p=orbitGeo["p"];
  e=orbitGeo["e"];
  x=orbitGeo["Inclination"];
  \[ScriptCapitalI]=ArcCos[x];
  Ehat=orbitGeo["Energy"];(* geodesic constants of motion *)
  Lzhat=orbitGeo["AngularMomentum"];
  Qhat=orbitGeo["CarterConstant"];
  Khat=Qhat+(Lzhat-a*Ehat)^2;
  \[CapitalUpsilon]rhat=orbitGeo["RadialFrequency"];(* geodesic frequencies *)
  \[CapitalUpsilon]\[Theta]hat=orbitGeo["PolarFrequency"];
  rhat[\[Lambda]_]:=rhat[\[Lambda]]=orbitGeo["Trajectory"][[2]][\[Lambda]];(* geodesic trajectory *)
  zhat[\[Lambda]_]:=zhat[\[Lambda]]=Cos[orbitGeo["Trajectory"][[3]][\[Lambda]]];
  \[CapitalLambda]rhat=2Pi/\[CapitalUpsilon]rhat;
  \[CapitalLambda]\[Theta]hat=2Pi/\[CapitalUpsilon]\[Theta]hat;
  stepsr=4*nmax;(* steps for numerical integration *)
  steps\[Theta]=4*kmax;
  ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];(* matrices of discrete Fourier transform *)
  ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/steps\[Theta]],Precision[{a,p,e,x}]],{j,1,steps\[Theta]},{k,-kmax,kmax}];
  Print["Calculating t, \[Phi] equations"];
  \[ScriptCapitalR]tlist=Table[\[ScriptCapitalR]tfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];(* list of values of \[ScriptCapitalR]t for Fourier transform *)
  \[ScriptCapitalR]\[Phi]list=Table[\[ScriptCapitalR]\[Phi]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalR]tfour=ExpniTable . \[ScriptCapitalR]tlist . ExpjkTable/(stepsr*steps\[Theta]);(* Discrete Fourier Transform of \[ScriptCapitalR]t and \[ScriptCapitalR]\[Phi] *)
  \[ScriptCapitalR]\[Phi]four=ExpniTable . \[ScriptCapitalR]\[Phi]list . ExpjkTable/(stepsr*steps\[Theta]);
  Print["Calculating r equation"];
  \[ScriptCapitalF]rlist=Table[\[ScriptCapitalF]rfunc[\[Chi]rhat[\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalF]\[GothicR]list=Table[\[ScriptCapitalF]\[GothicR]func[],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalG]rlist=Table[\[ScriptCapitalG]rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],\[Chi]rhat[\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalG]\[GothicR]list=Table[\[ScriptCapitalG]\[GothicR]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalG]\[Theta]list=Table[\[ScriptCapitalG]\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalG]\[GothicZ]list=Table[\[ScriptCapitalG]\[GothicZ]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalH]rlist=Table[\[ScriptCapitalH]rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]],d2rd\[Lambda]2[rhat[\[Lambda]r]],\[Chi]rhat[\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalH]\[GothicR]list=Table[\[ScriptCapitalH]\[GothicR]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalH]\[Theta]list=Table[\[ScriptCapitalH]\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalH]\[GothicZ]list=Table[\[ScriptCapitalH]\[GothicZ]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalI]1rlist=Table[\[ScriptCapitalI]1rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],d2rd\[Lambda]2[rhat[\[Lambda]r]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalI]1\[Theta]list=Table[\[ScriptCapitalI]1\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalI]2list=Table[\[ScriptCapitalI]2func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalI]3list=Table[\[ScriptCapitalI]3func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalJ]list=Table[\[ScriptCapitalJ]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],\[Delta]udst[\[Lambda]r,\[Lambda]\[Theta]],\[Delta]uds\[Phi][\[Lambda]r,\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalF]rfour=ExpniTable . \[ScriptCapitalF]rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalF]\[GothicR]four=ExpniTable . \[ScriptCapitalF]\[GothicR]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalG]rfour=ExpniTable . \[ScriptCapitalG]rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalG]\[GothicR]four=ExpniTable . \[ScriptCapitalG]\[GothicR]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalG]\[Theta]four=ExpniTable . \[ScriptCapitalG]\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalG]\[GothicZ]four=ExpniTable . \[ScriptCapitalG]\[GothicZ]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalH]rfour=ExpniTable . \[ScriptCapitalH]rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalH]\[GothicR]four=ExpniTable . \[ScriptCapitalH]\[GothicR]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalH]\[Theta]four=ExpniTable . \[ScriptCapitalH]\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalH]\[GothicZ]four=ExpniTable . \[ScriptCapitalH]\[GothicZ]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalI]1rfour=ExpniTable . \[ScriptCapitalI]1rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalI]1\[Theta]four=ExpniTable . \[ScriptCapitalI]1\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalI]2four=ExpniTable . \[ScriptCapitalI]2list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalI]3four=ExpniTable . \[ScriptCapitalI]3list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalJ]four=ExpniTable . \[ScriptCapitalJ]list . ExpjkTable/(stepsr*steps\[Theta]);
  Print["Calculating \[Theta] equation"];
  \[ScriptCapitalQ]\[Theta]list=Table[\[ScriptCapitalQ]\[Theta]func[zhat[\[Lambda]\[Theta]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalQ]\[GothicZ]list=Table[\[ScriptCapitalQ]\[GothicZ]func[zhat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalS]rlist=Table[\[ScriptCapitalS]rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],\[Chi]rhat[\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalS]\[GothicR]list=Table[\[ScriptCapitalS]\[GothicR]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalS]\[Theta]list=Table[\[ScriptCapitalS]\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalS]\[GothicZ]list=Table[\[ScriptCapitalS]\[GothicZ]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalT]rlist=Table[\[ScriptCapitalT]rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]],\[Chi]rhat[\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalT]\[GothicR]list=Table[\[ScriptCapitalT]\[GothicR]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalT]\[Theta]list=Table[\[ScriptCapitalT]\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]],d2zd\[Lambda]2[zhat[\[Lambda]\[Theta]]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalT]\[GothicZ]list=Table[\[ScriptCapitalT]\[GothicZ]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]],d2zd\[Lambda]2[zhat[\[Lambda]\[Theta]]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalU]1rlist=Table[\[ScriptCapitalU]1rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalU]1\[Theta]list=Table[\[ScriptCapitalU]1\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]],d2zd\[Lambda]2[zhat[\[Lambda]\[Theta]]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalU]2list=Table[\[ScriptCapitalU]2func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalU]3list=Table[\[ScriptCapitalU]3func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalV]list=Table[\[ScriptCapitalV]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],\[Delta]udst[\[Lambda]r,\[Lambda]\[Theta]],\[Delta]uds\[Phi][\[Lambda]r,\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalQ]\[Theta]four=ExpniTable . \[ScriptCapitalQ]\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalQ]\[GothicZ]four=ExpniTable . \[ScriptCapitalQ]\[GothicZ]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalS]rfour=ExpniTable . \[ScriptCapitalS]rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalS]\[GothicR]four=ExpniTable . \[ScriptCapitalS]\[GothicR]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalS]\[Theta]four=ExpniTable . \[ScriptCapitalS]\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalS]\[GothicZ]four=ExpniTable . \[ScriptCapitalS]\[GothicZ]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalT]rfour=ExpniTable . \[ScriptCapitalT]rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalT]\[GothicR]four=ExpniTable . \[ScriptCapitalT]\[GothicR]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalT]\[Theta]four=ExpniTable . \[ScriptCapitalT]\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalT]\[GothicZ]four=ExpniTable . \[ScriptCapitalT]\[GothicZ]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalU]1rfour=ExpniTable . \[ScriptCapitalU]1rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalU]1\[Theta]four=ExpniTable . \[ScriptCapitalU]1\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalU]2four=ExpniTable . \[ScriptCapitalU]2list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalU]3four=ExpniTable . \[ScriptCapitalU]3list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalV]four=ExpniTable . \[ScriptCapitalV]list . ExpjkTable/(stepsr*steps\[Theta]);
  Print["Calculating norm equation"];
  \[ScriptCapitalK]rlist=Table[\[ScriptCapitalK]rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],\[Chi]rhat[\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalK]\[GothicR]list=Table[\[ScriptCapitalK]\[GothicR]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalK]\[Theta]list=Table[\[ScriptCapitalK]\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalK]\[GothicZ]list=Table[\[ScriptCapitalK]\[GothicZ]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalM]rlist=Table[\[ScriptCapitalM]rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]],\[Chi]rhat[\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalM]\[GothicR]list=Table[\[ScriptCapitalM]\[GothicR]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalM]\[Theta]list=Table[\[ScriptCapitalM]\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]],\[Chi]\[Theta]hat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalM]\[GothicZ]list=Table[\[ScriptCapitalM]\[GothicZ]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalN]1rlist=Table[\[ScriptCapitalN]1rfunc[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],drd\[Lambda][\[Lambda]r]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalN]1\[Theta]list=Table[\[ScriptCapitalN]1\[Theta]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],dzd\[Lambda][\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalN]2list=Table[\[ScriptCapitalN]2func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalN]3list=Table[\[ScriptCapitalN]3func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalP]list=Table[\[ScriptCapitalP]func[rhat[\[Lambda]r],zhat[\[Lambda]\[Theta]],\[Delta]udst[\[Lambda]r,\[Lambda]\[Theta]],\[Delta]uds\[Phi][\[Lambda]r,\[Lambda]\[Theta]]],{\[Lambda]r,\[CapitalLambda]rhat/2/stepsr,\[CapitalLambda]rhat,\[CapitalLambda]rhat/stepsr},{\[Lambda]\[Theta],\[CapitalLambda]\[Theta]hat/2/steps\[Theta],\[CapitalLambda]\[Theta]hat,\[CapitalLambda]\[Theta]hat/steps\[Theta]}];
  \[ScriptCapitalK]rfour=ExpniTable . \[ScriptCapitalK]rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalK]\[GothicR]four=ExpniTable . \[ScriptCapitalK]\[GothicR]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalK]\[Theta]four=ExpniTable . \[ScriptCapitalK]\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalK]\[GothicZ]four=ExpniTable . \[ScriptCapitalK]\[GothicZ]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalM]rfour=ExpniTable . \[ScriptCapitalM]rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalM]\[GothicR]four=ExpniTable . \[ScriptCapitalM]\[GothicR]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalM]\[Theta]four=ExpniTable . \[ScriptCapitalM]\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalM]\[GothicZ]four=ExpniTable . \[ScriptCapitalM]\[GothicZ]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalN]1rfour=ExpniTable . \[ScriptCapitalN]1rlist . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalN]1\[Theta]four=ExpniTable . \[ScriptCapitalN]1\[Theta]list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalN]2four=ExpniTable . \[ScriptCapitalN]2list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalN]3four=ExpniTable . \[ScriptCapitalN]3list . ExpjkTable/(stepsr*steps\[Theta]);
  \[ScriptCapitalP]four=ExpniTable . \[ScriptCapitalP]list . ExpjkTable/(stepsr*steps\[Theta]);
  Print["Calculating submatrices"];(* Submatrices of M *)
  M1r=Flatten[Table[If[-nmax<=n-nn<=nmax,-nn^2*\[CapitalUpsilon]rhat^2*\[ScriptCapitalF]rfour[[n-nn+nmax+1,k+kmax+1]]-I*nn*\[CapitalUpsilon]rhat*\[ScriptCapitalG]rfour[[n-nn+nmax+1,k+kmax+1]]+\[ScriptCapitalH]rfour[[n-nn+nmax+1,k+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{nn,Join[Range[-nmax,-1],Range[1,nmax]]}],{{1,2},{3}}];(* r equation *)
  M1\[Theta]=Flatten[Table[If[-kmax<=k-kk<=kmax,-I*kk*\[CapitalUpsilon]\[Theta]hat*\[ScriptCapitalG]\[Theta]four[[n+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalH]\[Theta]four[[n+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}],{{1,2},{3}}];
  M1\[GothicR]=Flatten[Table[If[-nmax<=n-nn<=nmax&&-kmax<=k-kk<=kmax,(-I*nn*\[CapitalUpsilon]rhat-I*kk*\[CapitalUpsilon]\[Theta]hat)^2*\[ScriptCapitalF]\[GothicR]four[[n-nn+nmax+1,k-kk+kmax+1]]-I*(nn*\[CapitalUpsilon]rhat+kk*\[CapitalUpsilon]\[Theta]hat)*\[ScriptCapitalG]\[GothicR]four[[n-nn+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalH]\[GothicR]four[[n-nn+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{nn,-nmax,nmax},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}],{{1,2},{3,4}}];
  M1\[GothicZ]=Flatten[Table[If[-nmax<=n-nn<=nmax&&-kmax<=k-kk<=kmax,-I*(nn*\[CapitalUpsilon]rhat+kk*\[CapitalUpsilon]\[Theta]hat)*\[ScriptCapitalG]\[GothicZ]four[[n-nn+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalH]\[GothicZ]four[[n-nn+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{nn,Join[Range[-nmax,-1],Range[1,nmax]]},{kk,-kmax+1,kmax-1,2}],{{1,2},{3,4}}];
  M1\[CapitalUpsilon]r=Flatten[\[ScriptCapitalI]1rfour[[All,1;;-1;;2]],{{1,2}}];(* column vectors as submatrices; even k-modes *)
  M1\[CapitalUpsilon]\[Theta]=Flatten[\[ScriptCapitalI]1\[Theta]four[[All,1;;-1;;2]],{{1,2}}];
  M1t=Flatten[\[ScriptCapitalI]2four[[All,1;;-1;;2]],{{1,2}}];
  M1\[Phi]=Flatten[\[ScriptCapitalI]3four[[All,1;;-1;;2]],{{1,2}}];
  M1=Transpose[Join[Transpose[-Im[M1r]],Transpose[-Im[M1\[Theta]]],Transpose[Re[M1\[GothicR]]],Transpose[Re[M1\[GothicZ]]],Re[{M1\[CapitalUpsilon]r,M1\[CapitalUpsilon]\[Theta],M1t,M1\[Phi]}]]];
  M2r=Flatten[Table[If[-nmax<=n-nn<=nmax,-I*nn*\[CapitalUpsilon]rhat*\[ScriptCapitalS]rfour[[n-nn+nmax+1,k+kmax+1]]+\[ScriptCapitalT]rfour[[n-nn+nmax+1,k+kmax+1]],0],{n,-nmax,nmax},{k,-kmax+1,kmax-1,2},{nn,Join[Range[-nmax,-1],Range[1,nmax]]}],{{1,2},{3}}];(* \[Theta] equation *)
  M2\[Theta]=Flatten[Table[If[-kmax<=k-kk<=kmax,-kk^2*\[CapitalUpsilon]\[Theta]hat^2*\[ScriptCapitalQ]\[Theta]four[[n+nmax+1,k-kk+kmax+1]]-I*kk*\[CapitalUpsilon]\[Theta]hat*\[ScriptCapitalS]\[Theta]four[[n+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalT]\[Theta]four[[n+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax+1,kmax-1,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}],{{1,2},{3}}];
  M2\[GothicR]=Flatten[Table[If[-nmax<=n-nn<=nmax&&-kmax<=k-kk<=kmax,-I*(nn*\[CapitalUpsilon]rhat+kk*\[CapitalUpsilon]\[Theta]hat)*\[ScriptCapitalS]\[GothicR]four[[n-nn+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalT]\[GothicR]four[[n-nn+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax+1,kmax-1,2},{nn,-nmax,nmax},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}],{{1,2},{3,4}}];
  M2\[GothicZ]=Flatten[Table[If[-nmax<=n-nn<=nmax&&-kmax<=k-kk<=kmax,(-I*nn*\[CapitalUpsilon]rhat-I*kk*\[CapitalUpsilon]\[Theta]hat)^2*\[ScriptCapitalQ]\[GothicZ]four[[n-nn+nmax+1,k-kk+kmax+1]]-I*(nn*\[CapitalUpsilon]rhat+kk*\[CapitalUpsilon]\[Theta]hat)*\[ScriptCapitalS]\[GothicZ]four[[n-nn+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalT]\[GothicZ]four[[n-nn+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax+1,kmax-1,2},{nn,Join[Range[-nmax,-1],Range[1,nmax]]},{kk,-kmax+1,kmax-1,2}],{{1,2},{3,4}}];
  M2\[CapitalUpsilon]r=Flatten[\[ScriptCapitalU]1rfour[[All,2;;-2;;2]],{{1,2}}];(* odd k-modes *)
  M2\[CapitalUpsilon]\[Theta]=Flatten[\[ScriptCapitalU]1\[Theta]four[[All,2;;-2;;2]],{{1,2}}];
  M2t=Flatten[\[ScriptCapitalU]2four[[All,2;;-2;;2]],{{1,2}}];
  M2\[Phi]=Flatten[\[ScriptCapitalU]3four[[All,2;;-2;;2]],{{1,2}}];
  M2=Transpose[Join[Transpose[-Im@M2r],Transpose[-Im@M2\[Theta]],Transpose[Re@M2\[GothicR]],Transpose[Re@M2\[GothicZ]],Re@{M2\[CapitalUpsilon]r,M2\[CapitalUpsilon]\[Theta],M2t,M2\[Phi]}]];
  M3r=Flatten[Table[If[-nmax<=n-nn<=nmax,-I*nn*\[CapitalUpsilon]rhat*\[ScriptCapitalK]rfour[[n-nn+nmax+1,k+kmax+1]]+\[ScriptCapitalM]rfour[[n-nn+nmax+1,k+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{nn,Join[Range[-nmax,-1],Range[1,nmax]]}],{{1,2},{3}}];(* norm equation *)
  M3\[Theta]=Flatten[Table[If[-kmax<=k-kk<=kmax,-I*kk*\[CapitalUpsilon]\[Theta]hat*\[ScriptCapitalK]\[Theta]four[[n+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalM]\[Theta]four[[n+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}],{{1,2},{3}}];
  M3\[GothicR]=Flatten[Table[If[-nmax<=n-nn<=nmax&&-kmax<=k-kk<=kmax,-I*(nn*\[CapitalUpsilon]rhat+kk*\[CapitalUpsilon]\[Theta]hat)*\[ScriptCapitalK]\[GothicR]four[[n-nn+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalM]\[GothicR]four[[n-nn+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{nn,-nmax,nmax},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}],{{1,2},{3,4}}];
  M3\[GothicZ]=Flatten[Table[If[-nmax<=n-nn<=nmax&&-kmax<=k-kk<=kmax,-I*(nn*\[CapitalUpsilon]rhat+kk*\[CapitalUpsilon]\[Theta]hat)*\[ScriptCapitalK]\[GothicZ]four[[n-nn+nmax+1,k-kk+kmax+1]]+\[ScriptCapitalM]\[GothicZ]four[[n-nn+nmax+1,k-kk+kmax+1]],0],{n,-nmax,nmax},{k,-kmax,kmax,2},{nn,Join[Range[-nmax,-1],Range[1,nmax]]},{kk,-kmax+1,kmax-1,2}],{{1,2},{3,4}}];
  M3\[CapitalUpsilon]r=Flatten[\[ScriptCapitalN]1rfour[[All,1;;-1;;2]],{{1,2}}];
  M3\[CapitalUpsilon]\[Theta]=Flatten[\[ScriptCapitalN]1\[Theta]four[[All,1;;-1;;2]],{{1,2}}];
  M3t=Flatten[\[ScriptCapitalN]2four[[All,1;;-1;;2]],{{1,2}}];
  M3\[Phi]=Flatten[\[ScriptCapitalN]3four[[All,1;;-1;;2]],{{1,2}}];
  M3=Transpose[Join[Transpose[-Im@M3r],Transpose[-Im@M3\[Theta]],Transpose[Re@M3\[GothicR]],Transpose[Re@M3\[GothicZ]],Re@{M3\[CapitalUpsilon]r,M3\[CapitalUpsilon]\[Theta],M3t,M3\[Phi]}]];
  M=Join[M1,M2,M3];
  c=Re[-Join[Flatten[\[ScriptCapitalJ]four[[All,1;;-1;;2]]],Flatten[\[ScriptCapitalV]four[[All,2;;-2;;2]]],Flatten[\[ScriptCapitalP]four[[All,1;;-1;;2]]]]];(* source term, even, odd, even k-modes for r, \[Theta], norm equations *)
  Print["Calculating v"];
  v=LeastSquares[M,c];
  \[Delta]\[Chi]rSfourIm=v[[1;;2*nmax]];(* Imaginary part \[Delta]\[Chi]rS_n *)
  \[Delta]\[Chi]\[Theta]SfourIm=v[[2*nmax+1;;2*nmax+kmax]];(* Imaginary part of even k-modes of \[Delta]\[Chi]\[Theta]S_k (odd k-modes are zero) *)
  \[Delta]\[GothicR]four=Partition[v[[2*nmax+kmax+1;;2*nmax+kmax+kmax*(2nmax+1)]],kmax];(* all n-modes and even k-modes of \[Delta]\[GothicR]_nk  *)
  \[Delta]\[GothicZ]four=Partition[v[[2*nmax+kmax+kmax*(2nmax+1)+1;;2*nmax+kmax+kmax*(2nmax+1)+2*nmax*kmax]],kmax];(* all n-modes and odd k-modes of \[Delta]\[GothicZ]_nk  *)
  \[CapitalUpsilon]rS=Re[v[[-4]]];(* corrections to the frequencies *)
  \[CapitalUpsilon]\[Theta]S=Re[v[[-3]]];
  udt0S=Re[v[[-2]]];
  ud\[Phi]0S=Re[v[[-1]]];
  <|
    "\[ScriptCapitalR]tfourIm"->Im[\[ScriptCapitalR]tfour],
    "\[ScriptCapitalR]\[Phi]fourIm"->Im[\[ScriptCapitalR]\[Phi]four],
    "\[Delta]\[Chi]rSfourIm"->Re[\[Delta]\[Chi]rSfourIm],
    "\[Delta]\[Chi]\[Theta]SfourIm"->Re[\[Delta]\[Chi]\[Theta]SfourIm],
    "\[Delta]\[GothicR]four"->Re[\[Delta]\[GothicR]four],
    "\[Delta]\[GothicZ]four"->Re[\[Delta]\[GothicZ]four],
    "\[CapitalUpsilon]rS"->Re[\[CapitalUpsilon]rS],
    "\[CapitalUpsilon]\[Theta]S"->Re[\[CapitalUpsilon]\[Theta]S],
    "udt0S"->Re[udt0S],
    "ud\[Phi]0S"->Re[ud\[Phi]0S],
    "error"->Norm[M . v-c],(* error estimate from the residual *)
    "nmax"->nmax,
    "kmax"->kmax
  |>
]


KerrSpinOrbitCorrection2[a_,p_,e_,x_,correction_]:=Module[{\[ScriptCapitalI]=ArcCos[x],prec,\[ScriptCapitalR]\[Phi]fourIm,\[ScriptCapitalR]tfourIm,nmax,kmax,orbit,Ehat,Lzhat,Qhat,Khat,\[CapitalUpsilon]rhat,\[CapitalUpsilon]\[Theta]hat,\[CapitalUpsilon]\[Phi]hat,\[CapitalGamma]hat,
    \[CapitalDelta]tr,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]r,\[CapitalDelta]\[Phi]\[Theta],r1,r2,r3,r4,kr,z1,z2,kz,stepsr,steps\[Theta],\[Delta]udtS,\[Delta]ud\[Phi]S,\[Delta]\[Chi]rSfourIm,\[Delta]\[Chi]\[Theta]SfourIm,\[Delta]\[Chi]rS,\[Delta]\[Chi]\[Theta]S,d\[Delta]\[Chi]rS,d\[Delta]\[Chi]\[Theta]S,\[Delta]\[GothicR]four,\[Delta]\[GothicZ]four,\[Delta]\[GothicR],\[Delta]\[GothicZ],d\[Delta]\[GothicR],d\[Delta]\[GothicZ],\[CapitalUpsilon]rS,\[CapitalUpsilon]\[Theta]S,\[CapitalUpsilon]\[Phi]S,\[CapitalGamma]S,
    udt0S,ud\[Phi]0S,drd\[Lambda],dzd\[Lambda],d\[Theta]d\[Lambda],rhat,zhat,\[Chi]rhat,\[Chi]\[Theta]hat,St,Sr,Sz,S\[Phi],\[CapitalDelta],\[CapitalSigma],rS,zS,UtS,U\[Phi]S,UtSlist,U\[Phi]Slist,UtSfour,U\[Phi]Sfour,Im\[CapitalDelta]tSfour,Im\[CapitalDelta]\[Phi]Sfour,\[CapitalDelta]tS,\[CapitalDelta]\[Phi]S,Uthat,U\[Phi]hat,
    Uthatlist,U\[Phi]hatlist,Uthatfour,U\[Phi]hatfour,\[CapitalSigma]gtt,\[CapitalSigma]g\[Phi]\[Phi],\[CapitalSigma]gt\[Phi],d\[CapitalSigma]gttdr,d\[CapitalSigma]gttdz,d\[CapitalSigma]g\[Phi]\[Phi]dr,d\[CapitalSigma]g\[Phi]\[Phi]dz,d\[CapitalSigma]gt\[Phi]dr,d\[CapitalSigma]gt\[Phi]dz,Kkr,Kkz,h1,h2,h3,ES,LS,OrbitCorrection,
    ExpniTable,ExpjkTable,udtS,ud\[Phi]S},
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Khat (r^2-3 a^2 z^2)))/(Sqrt[Khat] (r^2+a^2 z^2)^3);(* functions for the four-velocity equations *)
  h2[r_,z_]:=1/(Sqrt[Khat] (r^2+a^2 z^2)^3) (-Ehat Lzhat r^6+a^4 Ehat Lzhat r^2 z^4+a^2 Ehat Lzhat r^4 (-2+z^2)-a^6 Ehat Lzhat z^4 (-2+z^2)+a^7 Ehat^2 z^4 (-1+z^2)+a r^3 (Lzhat^2 r+Khat (-1+z^2)+Khat r (-1+2 z^2)+r^3 (z^2-Ehat^2 (-1+z^2)))+a^3 (-Ehat^2 r^4 (-1+z^2)+r z^2 (2 r^3+2 Khat r z^2-3 Khat (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Khat-Lzhat^2+r ((-1+z^2)+r (2-z^2+Ehat^2 (-1+z^2)))));
  h3[r_,z_]:=-((2 a r z)/(Sqrt[Khat] (r^2+a^2 z^2)^2));
  \[Chi]rhat[qr_]:=If[Mod[qr,Pi]==0,0,Floor[qr/(2Pi)]*2Pi+If[Mod[qr,2Pi]<Pi,ArcCos[(p-r3)/(e r3)*((e r3/(p-r3)-(p-r3)/(e r3))/((p-r3)/(e r3)+(1-2 JacobiSN[Kkr/Pi*qr,kr]^2))+1)],2Pi-ArcCos[(p-r3)/(e r3)*((e r3/(p-r3)-(p-r3)/(e r3))/((p-r3)/(e r3)+(1-2 JacobiSN[Kkr/Pi*qr,kr]^2))+1)]]];
  \[Chi]\[Theta]hat[qz_]:=Floor[qz/(2Pi)]*2Pi+If[Mod[qz,2Pi]<Pi,ArcCos[JacobiSN[Kkz*2*(qz+Pi/2)/Pi,kz]],2Pi-ArcCos[JacobiSN[Kkz*2*(qz+Pi/2)/Pi,kz]]];
  \[Delta]\[Chi]rS[wr_]:=2*Sum[\[Delta]\[Chi]rSfourIm[[n+nmax]]*Sin[n*wr],{n,1,nmax}];
  d\[Delta]\[Chi]rS[wr_]:=2*Sum[n*\[CapitalUpsilon]rhat*\[Delta]\[Chi]rSfourIm[[n+nmax]]*Cos[n*wr],{n,1,nmax}];
  \[Delta]\[Chi]\[Theta]S[wz_]:=2*Sum[\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*Sin[2k*wz],{k,1,kmax/2}];
  d\[Delta]\[Chi]\[Theta]S[wz_]:=2*Sum[2*k*\[CapitalUpsilon]\[Theta]hat*\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*Cos[2k*wz],{k,1,kmax/2}];
  \[Delta]\[GothicR][wr_,wz_]:=2*Sum[\[Delta]\[GothicR]four[[n+nmax+1,k+kmax/2]]*Cos[n*wr+2k*wz],{n,-nmax,nmax},{k,1,kmax/2}];
  d\[Delta]\[GothicR][wr_,wz_]:=-2*Sum[(n*\[CapitalUpsilon]rhat+2k*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicR]four[[n+nmax+1,k+kmax/2]]*Sin[n*wr+2k*wz],{n,-nmax,nmax},{k,1,kmax/2}];
  \[Delta]\[GothicZ][wr_,wz_]:=2*Sum[\[Delta]\[GothicZ]four[[n+nmax,k+kmax/2+1]]*Cos[n*wr+(2k+1)*wz],{n,1,nmax},{k,-kmax/2,kmax/2-1}];
  d\[Delta]\[GothicZ][wr_,wz_]:=-2*Sum[(n*\[CapitalUpsilon]rhat+(2k+1)*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicZ]four[[n+nmax,k+kmax/2+1]]*Sin[n*wr+(2k+1)*wz],{n,1,nmax},{k,-kmax/2,kmax/2-1}];
  drd\[Lambda][wr_]:=Sign[Pi-Mod[wr,2Pi]]*Sqrt[(Ehat*(rhat[wr]^2+a^2)-a*Lzhat)^2-\[CapitalDelta][rhat[wr]]*(rhat[wr]^2+Khat)];
  dzd\[Lambda][wz_]:=-Sign[Pi-Mod[wz,2Pi]]*Sqrt[Qhat-zhat[wz]^2*(a^2*(1-Ehat^2)*(1-zhat[wz]^2)+Lzhat^2+Qhat)];
  d\[Theta]d\[Lambda][wz_]:=-dzd\[Lambda][wz]/Sqrt[1-zhat[wz]^2];
  rhat[wr_]:=orbit["Trajectory"][[2]][wr];
  zhat[wz_]:=Cos[orbit["Trajectory"][[3]][wz]];
  St[wr_,wz_]:=-1/Sqrt[Khat]*(-a*(rhat[wr]*dzd\[Lambda][wz]+zhat[wz]*drd\[Lambda][wr])/\[CapitalSigma][rhat[wr],zhat[wz]]);
  Sr[wr_,wz_]:=-1/Sqrt[Khat]*(a*zhat[wz]*((rhat[wr]^2+a^2)*Ehat-a*Lzhat)/\[CapitalDelta][rhat[wr]]);
  Sz[wr_,wz_]:=-1/Sqrt[Khat]*(a*rhat[wr]*Ehat-(rhat[wr]*Lzhat)/(1-zhat[wz]^2));
  S\[Phi][wr_,wz_]:=-1/Sqrt[Khat]*((a^2*zhat[wz]*(1-zhat[wz]^2)*drd\[Lambda][wr]+rhat[wr]*(rhat[wr]^2+a^2)*dzd\[Lambda][wz])/\[CapitalSigma][rhat[wr],zhat[wz]]);
  \[CapitalDelta][r_]:=r^2-2*r+a^2;
  \[CapitalSigma][r_,z_]:=r^2+a^2*z^2;
  \[CapitalSigma]gtt[r_,z_]:=-((r^2+a^2)^2/(a^2+(-2+r) r))+a^2 (1-z^2);
  \[CapitalSigma]g\[Phi]\[Phi][r_,z_]:=1/(1-z^2)-a^2/(a^2+r (-2+r));
  \[CapitalSigma]gt\[Phi][r_]:=a(1-(r^2+a^2)/(a^2-2 r+r^2));
  d\[CapitalSigma]gttdr[r_]:=2 ((r^2+a^2) ((r^2+a^2) (r-1)-2 (r^2-2r+a^2) r))/(r^2-2r+a^2)^2;
  d\[CapitalSigma]gttdz[z_]:=-2 a^2 z;
  d\[CapitalSigma]gt\[Phi]dr[r_]:=2a ((r^2+a^2) (r-1)-(r^2-2r+a^2) r)/(r^2-2r+a^2)^2;
  d\[CapitalSigma]g\[Phi]\[Phi]dr[r_]:=(2 a^2 (r-1))/(a^2+(-2+r) r)^2;
  d\[CapitalSigma]g\[Phi]\[Phi]dz[z_]:=(2 z)/(-1+z^2)^2;
  rS[wr_,wz_]:=e*p*\[Delta]\[Chi]rS[wr]*Sin[\[Chi]rhat[wr]]/(1+e*Cos[\[Chi]rhat[wr]])^2+\[Delta]\[GothicR][wr,wz];
  zS[wr_,wz_]:=-Sin[\[ScriptCapitalI]]*\[Delta]\[Chi]\[Theta]S[wz]*Sin[\[Chi]\[Theta]hat[wz]]+\[Delta]\[GothicZ][wr,wz];
  Uthat[wr_,wz_]:=\[CapitalSigma]gtt[rhat[wr],zhat[wz]]*(-Ehat )+\[CapitalSigma]gt\[Phi][rhat[wr]]*(Lzhat);
  U\[Phi]hat[wr_,wz_]:=\[CapitalSigma]g\[Phi]\[Phi][rhat[wr],zhat[wz]]*(Lzhat)+\[CapitalSigma]gt\[Phi][rhat[wr]]*(-Ehat);
  \[CapitalDelta]tS[wr_,wz_]:=Sum[2*Im\[CapitalDelta]tSfour[[nmax+1,k+kmax+1]]*Sin[k*wz],{k,1,kmax}]+Sum[2*Im\[CapitalDelta]tSfour[[n+nmax+1,k+kmax+1]]*Sin[n*wr+k*wz],{n,1,nmax},{k,-kmax,kmax}];
  \[CapitalDelta]\[Phi]S[wr_,wz_]:=Sum[2*Im\[CapitalDelta]\[Phi]Sfour[[nmax+1,k+kmax+1]]*Sin[k*wz],{k,1,kmax}]+Sum[2*Im\[CapitalDelta]\[Phi]Sfour[[n+nmax+1,k+kmax+1]]*Sin[n*wr+k*wz],{n,1,nmax},{k,-kmax,kmax}];
  OrbitCorrection[wr_,wz_]:=Module[{sinkwz,coskwz,sinnwr,cosnwr,\[Delta]\[Chi]rSp,\[Chi]rhatp,\[Delta]\[Chi]\[Theta]Sp,\[Chi]\[Theta]hatp,k,n},
    sinkwz={Sin[wz]};
    coskwz={Cos[wz]};
    For[k=2,k<=kmax,k++,
      AppendTo[sinkwz,sinkwz[[k-1]]*coskwz[[1]]+coskwz[[k-1]]*sinkwz[[1]]];
      AppendTo[coskwz,coskwz[[k-1]]*coskwz[[1]]-sinkwz[[k-1]]*sinkwz[[1]]];
    ];
    sinnwr={Sin[wr]};
    cosnwr={Cos[wr]};
    For[n=2,n<=nmax,n++,
      AppendTo[sinnwr,sinnwr[[n-1]]*cosnwr[[1]]+cosnwr[[n-1]]*sinnwr[[1]]];
      AppendTo[cosnwr,cosnwr[[n-1]]*cosnwr[[1]]-sinnwr[[n-1]]*sinnwr[[1]]];
    ];
    Clear[n,k];
    \[Delta]\[Chi]rSp=(2*Sum[\[Delta]\[Chi]rSfourIm[[n+nmax]]*sinnwr[[n]],{n,1,nmax}]);
    \[Chi]rhatp=\[Chi]rhat[wr];
    \[Delta]\[Chi]\[Theta]Sp=(2*Sum[\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*sinkwz[[2k]],{k,1,kmax/2}]);
    \[Chi]\[Theta]hatp=\[Chi]\[Theta]hat[wz];
    <|
      "\[CapitalDelta]tS"->2*(Sum[Im\[CapitalDelta]tSfour[[nmax+1,k+kmax+1]]*sinkwz[[k]],{k,1,kmax}]+Sum[Im\[CapitalDelta]tSfour[[n+nmax+1,kmax+1]]*sinnwr[[n]],{n,1,nmax}]+Sum[Im\[CapitalDelta]tSfour[[n+nmax+1,k+kmax+1]]*(coskwz[[k]] sinnwr[[n]]+cosnwr[[n]] sinkwz[[k]])+Im\[CapitalDelta]tSfour[[n+nmax+1,-k+kmax+1]]*(coskwz[[k]] sinnwr[[n]]-cosnwr[[n]] sinkwz[[k]]),{n,1,nmax},{k,1,kmax}]),
      "rS"->e*p*\[Delta]\[Chi]rSp*Sin[\[Chi]rhatp]/(1+e*Cos[\[Chi]rhatp])^2+2*(Sum[\[Delta]\[GothicR]four[[nmax+1,k+kmax/2]]*coskwz[[2k]],{k,1,kmax/2}]+Sum[\[Delta]\[GothicR]four[[n+nmax+1,k+kmax/2]]*(cosnwr[[n]] coskwz[[2k]]-sinnwr[[n]] sinkwz[[2k]])+\[Delta]\[GothicR]four[[-n+nmax+1,k+kmax/2]]*(cosnwr[[n]] coskwz[[2k]]+sinnwr[[n]] sinkwz[[2k]]),{n,1,nmax},{k,1,kmax/2}]),
      "zS"->-Sin[\[ScriptCapitalI]]*\[Delta]\[Chi]\[Theta]Sp*Sin[\[Chi]\[Theta]hatp]+2*(Sum[\[Delta]\[GothicZ]four[[n+nmax,k+kmax/2+1]]*(cosnwr[[n]] coskwz[[2k+1]]-sinnwr[[n]] sinkwz[[2k+1]]),{n,1,nmax},{k,0,kmax/2-1}]+Sum[\[Delta]\[GothicZ]four[[n+nmax,-k+kmax/2+1]]*(cosnwr[[n]] coskwz[[2k-1]]+sinnwr[[n]] sinkwz[[2k-1]]),{n,1,nmax},{k,1,kmax/2}]),
      "\[CapitalDelta]\[Phi]S"->2*(Sum[Im\[CapitalDelta]\[Phi]Sfour[[nmax+1,k+kmax+1]]*sinkwz[[k]],{k,1,kmax}]+Sum[Im\[CapitalDelta]\[Phi]Sfour[[n+nmax+1,kmax+1]]*sinnwr[[n]],{n,1,nmax}]+Sum[Im\[CapitalDelta]\[Phi]Sfour[[n+nmax+1,k+kmax+1]]*(coskwz[[k]] sinnwr[[n]]+cosnwr[[n]] sinkwz[[k]])+Im\[CapitalDelta]\[Phi]Sfour[[n+nmax+1,-k+kmax+1]]*(coskwz[[k]] sinnwr[[n]]-cosnwr[[n]] sinkwz[[k]]),{n,1,nmax},{k,1,kmax}]),
      "UrS"->(p*e*(\[Chi]rhat'[wr]*(Cos[\[Chi]rhatp]*\[Delta]\[Chi]rSp*\[CapitalUpsilon]rhat+Sin[\[Chi]rhatp]*\[CapitalUpsilon]rS+2*e*Sin[\[Chi]rhatp]^2*\[Delta]\[Chi]rSp*\[CapitalUpsilon]rhat/(1+e*Cos[\[Chi]rhatp]))+Sin[\[Chi]rhatp]*(2*Sum[n*\[CapitalUpsilon]rhat*\[Delta]\[Chi]rSfourIm[[n+nmax]]*cosnwr[[n]],{n,1,nmax}]))/(1+e*Cos[\[Chi]rhatp])^2-2*(Sum[(2k*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicR]four[[nmax+1,k+kmax/2]]*sinkwz[[2k]],{k,1,kmax/2}]+Sum[(n*\[CapitalUpsilon]rhat+2k*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicR]four[[n+nmax+1,k+kmax/2]]*(coskwz[[2k]] sinnwr[[n]]+cosnwr[[n]] sinkwz[[2k]])+(-n*\[CapitalUpsilon]rhat+2k*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicR]four[[-n+nmax+1,k+kmax/2]]*(-coskwz[[2k]] sinnwr[[n]]+cosnwr[[n]] sinkwz[[2k]]),{n,1,nmax},{k,1,kmax/2}])),
      "UzS"->(-Sin[\[ScriptCapitalI]]*(\[Chi]\[Theta]hat'[wz]*(Cos[\[Chi]\[Theta]hatp]*\[Delta]\[Chi]\[Theta]Sp*\[CapitalUpsilon]\[Theta]hat+Sin[\[Chi]\[Theta]hatp]*\[CapitalUpsilon]\[Theta]S)+Sin[\[Chi]\[Theta]hatp]*(2*Sum[2*k*\[CapitalUpsilon]\[Theta]hat*\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*coskwz[[2k]],{k,1,kmax/2}]))-2*(Sum[(n*\[CapitalUpsilon]rhat+(2k+1)*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicZ]four[[n+nmax,k+kmax/2+1]]*(coskwz[[2k+1]] sinnwr[[n]]+cosnwr[[n]] sinkwz[[2k+1]]),{n,1,nmax},{k,0,kmax/2-1}]+Sum[(n*\[CapitalUpsilon]rhat+(-2k+1)*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicZ]four[[n+nmax,-k+kmax/2+1]]*(coskwz[[2k-1]] sinnwr[[n]]-cosnwr[[n]] sinkwz[[2k-1]]),{n,1,nmax},{k,1,kmax/2}]))
    |>
  ];
  UtS[wr_,wz_]:=((d\[CapitalSigma]gttdr[rhat[wr]]*rS[wr,wz]+d\[CapitalSigma]gttdz[zhat[wz]]*zS[wr,wz])*(-Ehat)+\[CapitalSigma]gtt[rhat[wr],zhat[wz]]*(udtS[wr,wz])+d\[CapitalSigma]gt\[Phi]dr[rhat[wr]]*rS[wr,wz]*(Lzhat)+\[CapitalSigma]gt\[Phi][rhat[wr]]*(ud\[Phi]S[wr,wz]));
  U\[Phi]S[wr_,wz_]:=((d\[CapitalSigma]g\[Phi]\[Phi]dr[rhat[wr]]*rS[wr,wz]+d\[CapitalSigma]g\[Phi]\[Phi]dz[zhat[wz]]*zS[wr,wz])*(Lzhat)+\[CapitalSigma]g\[Phi]\[Phi][rhat[wr],zhat[wz]]*(ud\[Phi]S[wr,wz])+d\[CapitalSigma]gt\[Phi]dr[rhat[wr]]*rS[wr,wz]*(-Ehat)+\[CapitalSigma]gt\[Phi][rhat[wr]]*(udtS[wr,wz]));
(****************************************)
  prec=Precision[{a,p,e,x}];
  orbit=KerrGeodesics`KerrGeoOrbit`KerrGeoOrbit[a,p,e,x,"Parametrization"->"Phases","Method"->"Analytic"];
  {Ehat,Lzhat,Qhat}=Values[orbit["ConstantsOfMotion"]];
  Khat=Qhat+(a*Ehat-Lzhat)^2;
  \[CapitalUpsilon]rhat=orbit["RadialFrequency"];
  \[CapitalUpsilon]\[Theta]hat=orbit["PolarFrequency"];
  \[CapitalUpsilon]\[Phi]hat=orbit["AzimuthalFrequency"];
  \[CapitalGamma]hat=orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"];
  \[CapitalDelta]tr=orbit["TrajectoryDeltas"]["\[CapitalDelta]tr"];
  \[CapitalDelta]t\[Theta]=orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"];
  \[CapitalDelta]\[Phi]r=orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]r"];
  \[CapitalDelta]\[Phi]\[Theta]=orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"];
  {r1,r2,r3,r4}=orbit["RadialRoots"];
  kr=(r1-r2)*(r3-r4)/(r1-r3)/(r2-r4);
  Kkr=EllipticK[kr];
  z1=Sin[\[ScriptCapitalI]];
  z2=Sqrt[a^2*(1-Ehat^2)+Lzhat^2/(1-z1^2)];
  kz=a^2*(1-Ehat^2)*z1^2/z2^2;
  Kkz=EllipticK[kz];
  nmax=correction["nmax"];
  kmax=correction["kmax"];
  stepsr=4*nmax;
  steps\[Theta]=4*kmax;
  \[ScriptCapitalR]\[Phi]fourIm=correction["\[ScriptCapitalR]\[Phi]fourIm"];
  \[ScriptCapitalR]tfourIm=correction["\[ScriptCapitalR]tfourIm"];
  \[Delta]\[Chi]rSfourIm=correction["\[Delta]\[Chi]rSfourIm"];
  \[Delta]\[Chi]\[Theta]SfourIm=correction["\[Delta]\[Chi]\[Theta]SfourIm"];
  \[Delta]\[GothicR]four=correction["\[Delta]\[GothicR]four"];
  \[Delta]\[GothicZ]four=correction["\[Delta]\[GothicZ]four"];
  \[CapitalUpsilon]rS=correction["\[CapitalUpsilon]rS"];
  \[CapitalUpsilon]\[Theta]S=correction["\[CapitalUpsilon]\[Theta]S"];
  udt0S=correction["udt0S"];
  ud\[Phi]0S=correction["ud\[Phi]0S"];
  Print["Calculating ES,LS"];
  ES=-udt0S+Sum[h1[rhat[wr],zhat[w\[Theta]]],{wr,2Pi/stepsr/2,2Pi,2Pi/stepsr},{w\[Theta],2Pi/steps\[Theta]/2,2Pi,2Pi/steps\[Theta]}]/stepsr/steps\[Theta];
  LS=ud\[Phi]0S-Sum[h2[rhat[wr],zhat[w\[Theta]]]+h3[rhat[wr],zhat[w\[Theta]]]*drd\[Lambda][wr]*dzd\[Lambda][w\[Theta]],{wr,2Pi/stepsr/2,2Pi,2Pi/stepsr},{w\[Theta],2Pi/steps\[Theta]/2,2Pi,2Pi/steps\[Theta]}]/stepsr/steps\[Theta];
  udtS[wr_,wz_]:=-(ES)+h1[rhat[wr],zhat[wz]];
  ud\[Phi]S[wr_,wz_]:=(LS)+(h2[rhat[wr],zhat[wz]]+h3[rhat[wr],zhat[wz]]*drd\[Lambda][wr]*dzd\[Lambda][wz]);
  Print["Calculating UtSlist, U\[Phi]Slist"];
  UtSlist=Table[UtS[wr,wz],{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)},{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  U\[Phi]Slist=Table[U\[Phi]S[wr,wz],{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)},{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  Print["Calculating Uthatlist, U\[Phi]hatlist"];
  Uthatlist=Table[Uthat[wr,wz],{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)},{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  U\[Phi]hatlist=Table[U\[Phi]hat[wr,wz],{wr,2Pi/(stepsr)/2,2Pi,2Pi/(stepsr)},{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  Print["Calculating UtSfour, U\[Phi]Sfour"];
  ExpniTable=Table[N[Exp[2Pi*I*n*(i-1/2)/stepsr],Precision[{a,p,e,x}]],{n,-nmax,nmax},{i,1,stepsr}];
  ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/steps\[Theta]],Precision[{a,p,e,x}]],{j,1,steps\[Theta]},{k,-kmax,kmax}];
  UtSfour=ExpniTable . UtSlist . ExpjkTable/(stepsr*steps\[Theta]);
  U\[Phi]Sfour=ExpniTable . U\[Phi]Slist . ExpjkTable/(stepsr*steps\[Theta]);
  Print["Calculating Uthatfour, U\[Phi]hatfour"];
  Uthatfour=ExpniTable . Uthatlist . ExpjkTable/(stepsr*steps\[Theta]);
  U\[Phi]hatfour=ExpniTable . U\[Phi]hatlist . ExpjkTable/(stepsr*steps\[Theta]);
  Print["Calculating Im\[CapitalDelta]tSfour, Im\[CapitalDelta]\[Phi]Sfour"];
  Im\[CapitalDelta]tSfour=Table[If[n==0&&k==0,0,Re[UtSfour[[n+nmax+1,k+kmax+1]]]/(n*\[CapitalUpsilon]rhat+k*\[CapitalUpsilon]\[Theta]hat)-Re[Uthatfour[[n+nmax+1,k+kmax+1]]]/(n*\[CapitalUpsilon]rhat+k*\[CapitalUpsilon]\[Theta]hat)^2*(n*\[CapitalUpsilon]rS+k*\[CapitalUpsilon]\[Theta]S)],{n,-nmax,nmax},{k,-kmax,kmax}];
  Im\[CapitalDelta]\[Phi]Sfour=Table[If[n==0&&k==0,0,Re[U\[Phi]Sfour[[n+nmax+1,k+kmax+1]]]/(n*\[CapitalUpsilon]rhat+k*\[CapitalUpsilon]\[Theta]hat)-Re[U\[Phi]hatfour[[n+nmax+1,k+kmax+1]]]/(n*\[CapitalUpsilon]rhat+k*\[CapitalUpsilon]\[Theta]hat)^2*(n*\[CapitalUpsilon]rS+k*\[CapitalUpsilon]\[Theta]S)],{n,-nmax,nmax},{k,-kmax,kmax}];
  \[CapitalGamma]S=Re[UtSfour[[nmax+1,kmax+1]]];
  \[CapitalUpsilon]\[Phi]S=Re[U\[Phi]Sfour[[nmax+1,kmax+1]]];
  <|
    "a"->a,
    "p"->p,
    "e"->e,
    "\[ScriptCapitalI]"->\[ScriptCapitalI],
    "Ehat"->Ehat,
    "Lzhat"->Lzhat,
    "Khat"->Khat,
    "TrajectoryCorrection"->{Function[{wr,wz},\[CapitalDelta]tS[wr,wz]],Function[{wr,wz},rS[wr,wz]],Function[{wr,wz},zS[wr,wz]],Function[{wr,wz},\[CapitalDelta]\[Phi]S[wr,wz]]},
    "TrajectoryGeo"->orbit["Trajectory"],
    "MinoFrequenciesGeo"->{\[CapitalGamma]hat,\[CapitalUpsilon]rhat,\[CapitalUpsilon]\[Theta]hat,\[CapitalUpsilon]\[Phi]hat},
    "BLFrequenciesGeo"->{\[CapitalUpsilon]rhat/\[CapitalGamma]hat,\[CapitalUpsilon]\[Theta]hat/\[CapitalGamma]hat,\[CapitalUpsilon]\[Phi]hat/\[CapitalGamma]hat},
    "FourVelocityCorrection"->{Function[{wr,wz},(udtS[wr,wz])],Function[{wr,wz},(p*e*(\[Chi]rhat'[wr]*(Cos[\[Chi]rhat[wr]]*\[Delta]\[Chi]rS[wr]*\[CapitalUpsilon]rhat+Sin[\[Chi]rhat[wr]]*\[CapitalUpsilon]rS+2*e*Sin[\[Chi]rhat[wr]]^2*\[Delta]\[Chi]rS[wr]*\[CapitalUpsilon]rhat/(1+e*Cos[\[Chi]rhat[wr]]))+Sin[\[Chi]rhat[wr]]*d\[Delta]\[Chi]rS[wr])/(1+e*Cos[\[Chi]rhat[wr]])^2+d\[Delta]\[GothicR][wr,wz])],Function[{wr,wz},(-Sin[\[ScriptCapitalI]]*(\[Chi]\[Theta]hat'[wz]*(Cos[\[Chi]\[Theta]hat[wz]]*\[Delta]\[Chi]\[Theta]S[wz]*\[CapitalUpsilon]\[Theta]hat+Sin[\[Chi]\[Theta]hat[wz]]*\[CapitalUpsilon]\[Theta]S)+Sin[\[Chi]\[Theta]hat[wz]]*d\[Delta]\[Chi]\[Theta]S[wz])+d\[Delta]\[GothicZ][wr,wz])],Function[{wr,wz},(ud\[Phi]S[wr,wz])]},
    "Spin4vector"->{Function[{wr,wz},St[wr,wz]],Function[{wr,wz},Sr[wr,wz]],Function[{wr,wz},Sz[wr,wz]],Function[{wr,wz},S\[Phi][wr,wz]]},
    "TrajectoryDeltas"->orbit["TrajectoryDeltas"],
    "MinoFrequenciesCorrection"->{\[CapitalGamma]S,\[CapitalUpsilon]rS,\[CapitalUpsilon]\[Theta]S,\[CapitalUpsilon]\[Phi]S},
    "BLFrequenciesCorrection"->{\[CapitalUpsilon]rS/\[CapitalGamma]hat-\[CapitalUpsilon]rhat/\[CapitalGamma]hat^2*\[CapitalGamma]S,\[CapitalUpsilon]\[Theta]S/\[CapitalGamma]hat-\[CapitalUpsilon]\[Theta]hat/\[CapitalGamma]hat^2*\[CapitalGamma]S,\[CapitalUpsilon]\[Phi]S/\[CapitalGamma]hat-\[CapitalUpsilon]\[Phi]hat/\[CapitalGamma]hat^2*\[CapitalGamma]S},
    "ES"->ES,
    "LS"->LS,
    "\[Delta]\[Chi]rS"->\[Delta]\[Chi]rS,
    "\[Delta]\[Chi]\[Theta]S"->\[Delta]\[Chi]\[Theta]S,
    "\[Delta]r"->\[Delta]\[GothicR],
    "\[Delta]z"->\[Delta]\[GothicZ],
    "OrbitCorrection"->Function[{wr,wz},OrbitCorrection[wr,wz]],
    "error"->correction["error"]
  |>
]


KerrSpinOrbitCorrection[a_, p_, e_, x_, nmax_?(IntegerQ[#] && # > 0&), kmax_?EvenQ]:=Module[{orbit,correction},
  orbit = KerrGeodesics`KerrGeoOrbit`KerrGeoOrbit[a, p, e, x, "Parametrization"->"Mino", "Method"->"Analytic"];
  correction = Correction[orbit, nmax, kmax];
  KerrSpinOrbitCorrection2[a, p, e, x, correction]
]


KerrSpinOrbitNum[orbitCorrection_,spar_]:=Module[{\[CapitalGamma],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]},
  {\[CapitalGamma],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]}=orbitCorrection["MinoFrequenciesGeo"]+spar*orbitCorrection["MinoFrequenciesCorrection"];
  \[CapitalOmega]r=\[CapitalUpsilon]r/\[CapitalGamma];
  \[CapitalOmega]\[Theta]=\[CapitalUpsilon]\[Theta]/\[CapitalGamma];
  \[CapitalOmega]\[Phi]=\[CapitalUpsilon]\[Phi]/\[CapitalGamma];
  <|
    "a"->orbitCorrection["a"],
    "p"->orbitCorrection["p"],
    "e"->orbitCorrection["e"],
    "\[ScriptCapitalI]"->orbitCorrection["\[ScriptCapitalI]"],
    "s"->spar,
    "Ehat"->orbitCorrection["Ehat"],
    "Lzhat"->orbitCorrection["Lzhat"],
    "Khat"->orbitCorrection["Khat"],
    "BLFrequencies"->{\[CapitalOmega]r,\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]},
    "TrajectoryGeo"->orbitCorrection["TrajectoryGeo"],
    "TrajectoryDeltas"->orbitCorrection["TrajectoryDeltas"],
    "MinoFrequencies"->{\[CapitalGamma],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]},
    "MinoFrequenciesGeo"->orbitCorrection["MinoFrequenciesGeo"],
    "MinoFrequenciesCorrection"->orbitCorrection["MinoFrequenciesCorrection"],
    "ES"->orbitCorrection["ES"],
    "LS"->orbitCorrection["LS"],
    "OrbitCorrection"->orbitCorrection["OrbitCorrection"]
  |>
]


(* ::Subsection::Closed:: *)
(*Numerical spherical orbit*)


CorrectionSpherical[orbitGeo_,kmax_?EvenQ]:=Module[{a,p,e,x,\[ScriptCapitalI],\[ScriptCapitalR]tfunc,\[ScriptCapitalR]\[Phi]func,\[ScriptCapitalG]\[Theta]func,\[ScriptCapitalH]\[GothicR]func,\[ScriptCapitalH]\[Theta]func,\[ScriptCapitalI]1\[Theta]func,\[ScriptCapitalI]2func,
    \[ScriptCapitalI]3func,\[ScriptCapitalJ]func,\[ScriptCapitalQ]\[Theta]func,\[ScriptCapitalS]\[Theta]func,\[ScriptCapitalT]\[Theta]func,\[ScriptCapitalU]1\[Theta]func,\[ScriptCapitalU]2func,\[ScriptCapitalU]3func,\[ScriptCapitalV]func,\[ScriptCapitalK]\[Theta]func,
    \[ScriptCapitalM]\[Theta]func,\[ScriptCapitalN]1\[Theta]func,\[ScriptCapitalN]2func,\[ScriptCapitalN]3func,\[ScriptCapitalP]func,Ehat,Lzhat,Qhat,Khat,\[CapitalUpsilon]\[Theta]hat,zhat,dzd\[Lambda],\[CapitalDelta],\[CapitalSigma],d2zd\[Lambda]2,Z,\[Chi]\[Theta]hat,\[Delta]udst,\[Delta]uds\[Phi],
    steps\[Theta],ExpjkTable,\[ScriptCapitalR]tlist,\[ScriptCapitalR]tfour,\[ScriptCapitalR]\[Phi]list,\[ScriptCapitalR]\[Phi]four,\[ScriptCapitalG]\[Theta]four,\[ScriptCapitalG]\[Theta]list,
    \[ScriptCapitalH]\[GothicR]four,\[ScriptCapitalH]\[GothicR]list,\[ScriptCapitalH]\[Theta]four,\[ScriptCapitalH]\[Theta]list,\[ScriptCapitalI]1\[Theta]four,\[ScriptCapitalI]1\[Theta]list,\[ScriptCapitalI]2four,\[ScriptCapitalI]2list,\[ScriptCapitalI]3four,\[ScriptCapitalI]3list,\[ScriptCapitalJ]four,\[ScriptCapitalJ]list,\[ScriptCapitalQ]\[Theta]four,\[ScriptCapitalQ]\[Theta]list,
    \[ScriptCapitalS]\[Theta]four,\[ScriptCapitalS]\[Theta]list,\[ScriptCapitalT]\[Theta]four,\[ScriptCapitalT]\[Theta]list,\[ScriptCapitalU]1\[Theta]four,\[ScriptCapitalU]1\[Theta]list,\[ScriptCapitalU]2four,\[ScriptCapitalU]2list,\[ScriptCapitalU]3four,
    \[ScriptCapitalU]3list,\[ScriptCapitalV]four,\[ScriptCapitalV]list,\[ScriptCapitalK]\[Theta]four,\[ScriptCapitalK]\[Theta]list,\[ScriptCapitalM]\[Theta]four,\[ScriptCapitalM]\[Theta]list,\[ScriptCapitalN]1\[Theta]four,
    \[ScriptCapitalN]1\[Theta]list,\[ScriptCapitalN]2four,\[ScriptCapitalN]2list,\[ScriptCapitalN]3four,\[ScriptCapitalN]3list,\[ScriptCapitalP]four,\[ScriptCapitalP]list,M1\[Theta],M1\[GothicR],M1\[CapitalUpsilon]\[Theta],M1t,M1\[Phi],M1,M2\[Theta],M2\[GothicR],M2\[CapitalUpsilon]\[Theta],M2t,M2\[Phi],M2,M3\[Theta],M3\[GothicR],M3\[CapitalUpsilon]\[Theta],M3t,M3\[Phi],M3,
    M,c,v,\[Delta]\[Chi]\[Theta]SfourIm,\[Delta]\[GothicR]four,\[CapitalUpsilon]\[Theta]S,udt0S,ud\[Phi]0S},
  dzd\[Lambda][qz_]:=-Sign[Pi-Mod[qz,2Pi]]*Sqrt[Qhat-zhat[qz]^2*(a^2*(1-Ehat^2)*(1-zhat[qz]^2)+Lzhat^2+Qhat)];(* geodesic polar velocity *)
  \[CapitalDelta][r_]:=r^2-2*r+a^2;
  \[CapitalSigma][r_,z_]:=r^2+a^2*z^2;
  d2zd\[Lambda]2[z_]:=(2 a^2 (1-Ehat^2) z^3-2 z (Lzhat^2+Qhat+a^2 (1-Ehat^2) (1-z^2)))/2;(* second derivative of z=cos(\[Theta]) wrt \[Lambda] *)
  Z[z_]:=Qhat-z^2*(a^2*(1-Ehat^2)*(1-z^2)+Lzhat^2+Qhat);
  \[Chi]\[Theta]hat[qz_]:=\[Chi]\[Theta]hat[qz]=If[Mod[qz,2Pi]<Pi,ArcCos[Cos[orbitGeo["Trajectory"][[3]][qz]]/Sin[\[ScriptCapitalI]]]+2Pi*Floor[qz/(2Pi)],2Pi+2Pi*Floor[qz/(2Pi)]-ArcCos[Cos[orbitGeo["Trajectory"][[3]][qz]]/Sin[\[ScriptCapitalI]]]];
  \[Delta]udst[qz_]:=\[Delta]udst[qz]=Re[Sum[If[k==0,0,Im[\[ScriptCapitalR]tfour[[k+kmax+1]]]/(-k*\[CapitalUpsilon]\[Theta]hat)*Exp[-I*k*qz]],{k,-kmax,kmax,2}]];(* definition of \[Delta]udst summed over even k for aligned spin *)
  \[Delta]uds\[Phi][qz_]:=\[Delta]uds\[Phi][qz]=Re[Sum[If[k==0,0,Im[\[ScriptCapitalR]\[Phi]four[[k+kmax+1]]]/(-k*\[CapitalUpsilon]\[Theta]hat)*Exp[-I*k*qz]],{k,-kmax,kmax,2}]];
  \[ScriptCapitalR]tfunc[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=-((2 a^2 dzd\[Lambda] r z (3 r^2 (2 Khat+r^2)-2 a^2 (3 Khat+4 r^2) z^2+a^4 z^4)+drd\[Lambda] (3 Khat r^4-9 a^2 r^2 (2 Khat+r^2) z^2+a^4 (3 Khat+14 r^2) z^4-a^6 z^6))/(Sqrt[Khat] (r^2+a^2 z^2)^4));(* functionf for t and \[Phi] equations *)
  \[ScriptCapitalR]\[Phi]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=-((a (-2 dzd\[Lambda] r (a^2+r^2) z (3 r^2 (2 Khat+r^2)-2 a^2 (3 Khat+4 r^2) z^2+a^4 z^4)+drd\[Lambda] (-1+z^2) (3 Khat r^4-9 a^2 r^2 (2 Khat+r^2) z^2+a^4 (3 Khat+14 r^2) z^4-a^6 z^6)))/(Sqrt[Khat] (r^2+a^2 z^2)^4));
  \[ScriptCapitalG]\[Theta]func[r_,z_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=-((2 dzd\[Lambda] r (a^2+(-2+r) r) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/((-1+z^2) (r^2+a^2 z^2)));
  \[ScriptCapitalH]\[GothicR]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_]:=1/2 (1/(a^2+(-2+r) r)^2 (2 a^2 (-2 a^2 Ehat^2+4 a Ehat Lzhat+(-2+a^2) Lzhat^2)+2 r (2 a^2 (-a^2 Ehat^2-4 a Ehat Lzhat+Lzhat^2)+r (a (16 a Ehat^2+4 Ehat Lzhat-a Lzhat^2)-2 Ehat^2 r (8+2 a^2-5 r+r^2))))+(Lzhat^2 (a^4-4 a^2 r+2 a^2 r^2+r^4))/((a^2+r^2)^2 (-1+z))-(Lzhat^2 (a^4-4 a^2 r+2 a^2 r^2+r^4))/((a^2+r^2)^2 (1+z))-(8 r^3 (a (a Ehat-Lzhat)+Ehat r^2)^2)/((a^2+r^2) (r^2+a^2 z^2)^2)+(8 r (a (a Ehat-Lzhat)+Ehat r^2) (a^3 (a Ehat-Lzhat)+3 a^2 Ehat r^2+2 Ehat r^4))/((a^2+r^2)^2 (r^2+a^2 z^2)))+drd\[Lambda]^2 ((2-a^2+(-2+r) r)/(a^2+(-2+r) r)^2+(2 r^2)/(r^2+a^2 z^2)^2-1/(r^2+a^2 z^2))+(dzd\[Lambda]^2 (-a^2 r^2+r^4+a^2 (a^2+r (-4+3 r)) z^2))/((-1+z^2) (r^2+a^2 z^2)^2);
  \[ScriptCapitalH]\[Theta]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=1/((-1+z^2)^2 (r^2+a^2 z^2)^2) 2 r z (Lzhat^2 r^4+2 a^2 Lzhat^2 r^2 z^2+a^4 Lzhat^2 z^4+2 r^3 (-a Ehat-Lzhat+a Ehat z^2) (-a Ehat+Lzhat+a Ehat z^2)+2 a^2 r (-a Ehat+Lzhat+a Ehat z^2) (-a Ehat+Lzhat+(a Ehat-2 Lzhat) z^2)) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat]-(2 a^2 drd\[Lambda]^2 r z Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/(r^2+a^2 z^2)^2-(2 dzd\[Lambda]^2 r (a^2+(-2+r) r) z (r^2 Cot[\[ScriptCapitalI]]^2+a^2 (1-2 z^2+z^4 Csc[\[ScriptCapitalI]]^2)) Csc[\[Chi]\[Theta]hat] Sin[\[ScriptCapitalI]])/((-1+z^2)^2 (r^2+a^2 z^2)^2);
  \[ScriptCapitalI]1\[Theta]func[r_,z_,dzd\[Lambda]_]:=(2 dzd\[Lambda]^2 r (a^2+(-2+r) r)/\[CapitalUpsilon]\[Theta]hat)/((-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalI]2func[r_,z_]:=1/(a^2+(-2+r) r) 2 (a^3 (a Ehat-Lzhat)+r^2 (a (2 a Ehat+Lzhat)+Ehat (-4+r) r))-(4 r^2 (a (a Ehat-Lzhat)+Ehat r^2))/(r^2+a^2 z^2);
  \[ScriptCapitalI]3func[r_,z_]:=(2 (Lzhat (-2+r)^2 r^3-a Ehat r^3 (-4+3 r) (-1+z^2)+a^5 Ehat z^2 (-1+z^2)+a^4 Lzhat z^2 (1+(-1+r) z^2)+a^2 Lzhat r^2 (-1+(-3+2 r) z^2)-a^3 Ehat r^2 (-1+z^4)))/((a^2+(-2+r) r) (-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalJ]func[r_,z_,\[Delta]udst_,\[Delta]uds\[Phi]_]:=((a^2 Ehat-a Lzhat+Ehat r^2) (-3 Khat r^4+9 a^2 r^2 (2 Khat+r^2) z^2-a^4 (3 Khat+14 r^2) z^4+a^6 z^6))/(Sqrt[Khat] (r^2+a^2 z^2)^3)+\[ScriptCapitalI]2func[r,z]*\[Delta]udst+\[ScriptCapitalI]3func[r,z]*\[Delta]uds\[Phi];
  \[ScriptCapitalK]\[Theta]func[r_,z_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=(2 dzd\[Lambda] Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/((-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalM]\[Theta]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=(2 a^2 drd\[Lambda]^2 z Sin[\[Chi]\[Theta]hat] Sin[\[ScriptCapitalI]])/((a^2+(-2+r) r) (r^2+a^2 z^2)^2)-1/((a^2+(-2+r) r) (-1+z^2)^2 (r^2+a^2 z^2)^2) 2 z (Lzhat^2 r^4+2 a^2 Lzhat^2 r^2 z^2+a^4 Lzhat^2 z^4+2 r^3 (-a Ehat-Lzhat+a Ehat z^2) (-a Ehat+Lzhat+a Ehat z^2)+2 a^2 r (-a Ehat+Lzhat+a Ehat z^2) (-a Ehat+Lzhat+(a Ehat-2 Lzhat) z^2)) Sin[\[Chi]\[Theta]hat] Sin[\[ScriptCapitalI]]+(2 dzd\[Lambda]^2 z (r^2 Cot[\[ScriptCapitalI]]^2+a^2 (1-2 z^2+z^4 Csc[\[ScriptCapitalI]]^2)) Sin[\[ScriptCapitalI]])/((-1+z^2)^2 (r^2+a^2 z^2)^2 Sin[\[Chi]\[Theta]hat]);
  \[ScriptCapitalN]1\[Theta]func[r_,z_,dzd\[Lambda]_]:=-((2 dzd\[Lambda]^2/\[CapitalUpsilon]\[Theta]hat)/((-1+z^2) (r^2+a^2 z^2)));
  \[ScriptCapitalN]2func[r_,z_]:=2 Ehat+(4 r (a (a Ehat-Lzhat)+Ehat r^2))/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalN]3func[r_,z_]:=(-2 r (2 a Ehat-2 Lzhat+Lzhat r)-2 a (a Lzhat-2 Ehat r) z^2)/((a^2+(-2+r) r) (-1+z^2) (r^2+a^2 z^2));
  \[ScriptCapitalP]func[r_,z_,\[Delta]udst_,\[Delta]uds\[Phi]_]:=\[ScriptCapitalN]2func[r,z]*\[Delta]udst+\[ScriptCapitalN]3func[r,z]*\[Delta]uds\[Phi];
  \[ScriptCapitalQ]\[Theta]func[z_,\[Chi]\[Theta]hat_]:=(Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/Sqrt[1-z^2];(* functions for \[Theta] equation *)
  \[ScriptCapitalS]\[Theta]func[r_,z_,dzd\[Lambda]_,\[Chi]\[Theta]hat_]:=-((2 z Sqrt[1-z^2] (r^2 Cot[\[ScriptCapitalI]]^2+a^2 (1-2 z^2+z^4 Csc[\[ScriptCapitalI]]^2)) dzd\[Lambda] Sin[\[ScriptCapitalI]])/((-1+z^2)^2 (r^2+a^2 z^2) Sin[\[Chi]\[Theta]hat]));
  \[ScriptCapitalT]\[Theta]func[r_,z_,drd\[Lambda]_,dzd\[Lambda]_,d2zd\[Lambda]2_,\[Chi]\[Theta]hat_]:=1/((a^2+(-2+r) r) (1-z^2)^(5/2) (r^2+a^2 z^2)^2) (3 a^6 Lzhat^2 z^6+Lzhat^2 r^6 (1+2 z^2)+a^4 Lzhat^2 r^2 z^4 (7+2 z^2)+a^2 Lzhat^2 r^4 z^2 (5+4 z^2)+2 a^4 r z^2 (-(-a Ehat+Lzhat)^2+2 (-a Ehat+Lzhat)^2 z^2-(a^2 Ehat^2-2 a Ehat Lzhat+4 Lzhat^2) z^4)+2 a^2 r^3 ((-a Ehat+Lzhat)^2-(5 a Ehat-3 Lzhat) (a Ehat-Lzhat) z^2+(7 a^2 Ehat^2-10 a Ehat Lzhat-4 Lzhat^2) z^4+a Ehat (-3 a Ehat+4 Lzhat) z^6)-2 r^5 (-a^2 Ehat^2+Lzhat^2+2 (2 a^2 Ehat^2+Lzhat^2) z^2+a^2 Ehat^2 z^4 (-5+2 z^2))) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat]+(a^2 drd\[Lambda]^2 (a^2 z^2+r^2 (-1+2 z^2)) Sin[\[ScriptCapitalI]] Sin[\[Chi]\[Theta]hat])/((a^2+(-2+r) r) Sqrt[1-z^2] (r^2+a^2 z^2)^2)-1/(2 (1-z^2)^(5/2) (r^2+a^2 z^2)^2) dzd\[Lambda]^2 Csc[\[ScriptCapitalI]] Csc[\[Chi]\[Theta]hat] (2 r^4 (1+2 z^2) Cos[\[ScriptCapitalI]]^2+a^2 r^2 (-4 z^6+2 z^4 (7+2 Cos[2 \[ScriptCapitalI]])+z^2 (-5+3 Cos[2 \[ScriptCapitalI]])+2 Sin[\[ScriptCapitalI]]^2)+2 a^4 z^2 (z^4 (1+2 Cos[2 \[ScriptCapitalI]])-Sin[\[ScriptCapitalI]]^2+2 z^2 Sin[\[ScriptCapitalI]]^2))+(2 z Cos[\[ScriptCapitalI]] Cot[\[ScriptCapitalI]] Csc[\[Chi]\[Theta]hat] (dzd\[Lambda]^2 z+d2zd\[Lambda]2 (-z^2+Sin[\[ScriptCapitalI]]^2)))/((1-z^2)^(3/2) (-1+2 z^2+Cos[2 \[ScriptCapitalI]]));
  \[ScriptCapitalU]1\[Theta]func[r_,z_,dzd\[Lambda]_,d2zd\[Lambda]2_,\[Chi]\[Theta]hat_]:=-((-2 d2zd\[Lambda]2 (-1+z^2) (r^2+a^2 z^2)+2 dzd\[Lambda]^2 z (r^2+a^2 (-1+2 z^2)))/((1-z^2)^(3/2) (r^2+a^2 z^2) \[CapitalUpsilon]\[Theta]hat));
  \[ScriptCapitalU]2func[r_,z_]:=(4 a^2 r (a^2 Ehat-a Lzhat+Ehat r^2) z Sqrt[1-z^2])/((a^2+(-2+r) r) (r^2+a^2 z^2));
  \[ScriptCapitalU]3func[r_,z_]:=-((2 z (Lzhat (-2+r) r^3+a^4 Lzhat z^4-2 a^3 Ehat r (-1+z^2)^2+2 a^2 Lzhat r (1+(-2+r) z^2)))/((a^2+(-2+r) r) (1-z^2)^(3/2) (r^2+a^2 z^2)));
  \[ScriptCapitalV]func[r_,z_,\[Delta]udst_,\[Delta]uds\[Phi]_]:=-((2 a r z (3 r^2 (2 Khat+r^2)-2 a^2 (3 Khat+4 r^2) z^2+a^4 z^4) (Lzhat+a Ehat (-1+z^2)))/(Sqrt[Khat] (r^2+a^2 z^2)^3))/Sqrt[1-z^2]+\[ScriptCapitalU]2func[r,z]*\[Delta]udst+\[ScriptCapitalU]3func[r,z]*\[Delta]uds\[Phi];
  a=orbitGeo["a"];(*orbital parameters*)
  p=orbitGeo["p"];
  e=orbitGeo["e"];
  x=orbitGeo["Inclination"];
  \[ScriptCapitalI]=ArcCos[x];
  Ehat=orbitGeo["Energy"];(* geodesic constants of motion *)
  Lzhat=orbitGeo["AngularMomentum"];
  Qhat=orbitGeo["CarterConstant"];
  Khat=Qhat+(Lzhat-a*Ehat)^2;
  \[CapitalUpsilon]\[Theta]hat=orbitGeo["PolarFrequency"];
  zhat[qz_]:=zhat[qz]=Cos[orbitGeo["Trajectory"][[3]][qz]];
  steps\[Theta]=4*kmax;
  ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/steps\[Theta]],Precision[{a,p,e,x}]],{j,1,steps\[Theta]},{k,-kmax,kmax}];
  Print["Calculating t, \[Phi] equations"];
  \[ScriptCapitalR]tlist=Table[\[ScriptCapitalR]tfunc[p,zhat[qz],0,dzd\[Lambda][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];(* list of values of \[ScriptCapitalR]t for Fourier transform *)
  \[ScriptCapitalR]\[Phi]list=Table[\[ScriptCapitalR]\[Phi]func[p,zhat[qz],0,dzd\[Lambda][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalR]tfour=\[ScriptCapitalR]tlist . ExpjkTable/(steps\[Theta]);(* Discrete Fourier Transform of \[ScriptCapitalR]t and \[ScriptCapitalR]\[Phi] *)
  \[ScriptCapitalR]\[Phi]four=\[ScriptCapitalR]\[Phi]list . ExpjkTable/(steps\[Theta]);
  Print["Calculating r equation"];
  \[ScriptCapitalG]\[Theta]list=Table[\[ScriptCapitalG]\[Theta]func[p,zhat[qz],dzd\[Lambda][qz],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalH]\[GothicR]list=Table[\[ScriptCapitalH]\[GothicR]func[p,zhat[qz],0,dzd\[Lambda][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalH]\[Theta]list=Table[\[ScriptCapitalH]\[Theta]func[p,zhat[qz],0,dzd\[Lambda][qz],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalI]1\[Theta]list=Table[\[ScriptCapitalI]1\[Theta]func[p,zhat[qz],dzd\[Lambda][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalI]2list=Table[\[ScriptCapitalI]2func[p,zhat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalI]3list=Table[\[ScriptCapitalI]3func[p,zhat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalJ]list=Table[\[ScriptCapitalJ]func[p,zhat[qz],\[Delta]udst[qz],\[Delta]uds\[Phi][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalG]\[Theta]four=\[ScriptCapitalG]\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalH]\[GothicR]four=\[ScriptCapitalH]\[GothicR]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalH]\[Theta]four=\[ScriptCapitalH]\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalI]1\[Theta]four=\[ScriptCapitalI]1\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalI]2four=\[ScriptCapitalI]2list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalI]3four=\[ScriptCapitalI]3list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalJ]four=\[ScriptCapitalJ]list . ExpjkTable/(steps\[Theta]);
  Print["Calculating \[Theta] equation"];
  \[ScriptCapitalQ]\[Theta]list=Table[\[ScriptCapitalQ]\[Theta]func[zhat[qz],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalS]\[Theta]list=Table[\[ScriptCapitalS]\[Theta]func[p,zhat[qz],dzd\[Lambda][qz],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalT]\[Theta]list=Table[\[ScriptCapitalT]\[Theta]func[p,zhat[qz],0,dzd\[Lambda][qz],d2zd\[Lambda]2[zhat[qz]],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalU]1\[Theta]list=Table[\[ScriptCapitalU]1\[Theta]func[p,zhat[qz],dzd\[Lambda][qz],d2zd\[Lambda]2[zhat[qz]],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalU]2list=Table[\[ScriptCapitalU]2func[p,zhat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalU]3list=Table[\[ScriptCapitalU]3func[p,zhat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalV]list=Table[\[ScriptCapitalV]func[p,zhat[qz],\[Delta]udst[qz],\[Delta]uds\[Phi][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalQ]\[Theta]four=\[ScriptCapitalQ]\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalS]\[Theta]four=\[ScriptCapitalS]\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalT]\[Theta]four=\[ScriptCapitalT]\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalU]1\[Theta]four=\[ScriptCapitalU]1\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalU]2four=\[ScriptCapitalU]2list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalU]3four=\[ScriptCapitalU]3list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalV]four=\[ScriptCapitalV]list . ExpjkTable/(steps\[Theta]);
  Print["Calculating norm equation"];
  \[ScriptCapitalK]\[Theta]list=Table[\[ScriptCapitalK]\[Theta]func[p,zhat[qz],dzd\[Lambda][qz],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalM]\[Theta]list=Table[\[ScriptCapitalM]\[Theta]func[p,zhat[qz],0,dzd\[Lambda][qz],\[Chi]\[Theta]hat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalN]1\[Theta]list=Table[\[ScriptCapitalN]1\[Theta]func[p,zhat[qz],dzd\[Lambda][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalN]2list=Table[\[ScriptCapitalN]2func[p,zhat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalN]3list=Table[\[ScriptCapitalN]3func[p,zhat[qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalP]list=Table[\[ScriptCapitalP]func[p,zhat[qz],\[Delta]udst[qz],\[Delta]uds\[Phi][qz]],{qz,2Pi/2/steps\[Theta],2Pi,2Pi/steps\[Theta]}];
  \[ScriptCapitalK]\[Theta]four=\[ScriptCapitalK]\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalM]\[Theta]four=\[ScriptCapitalM]\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalN]1\[Theta]four=\[ScriptCapitalN]1\[Theta]list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalN]2four=\[ScriptCapitalN]2list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalN]3four=\[ScriptCapitalN]3list . ExpjkTable/(steps\[Theta]);
  \[ScriptCapitalP]four=\[ScriptCapitalP]list . ExpjkTable/(steps\[Theta]);
  Print["Calculating submatrices"];(* Submatrices of M *)
  M1\[Theta]=Table[If[-kmax<=k-kk<=kmax,-I*kk*\[CapitalUpsilon]\[Theta]hat*\[ScriptCapitalG]\[Theta]four[[k-kk+kmax+1]]+\[ScriptCapitalH]\[Theta]four[[k-kk+kmax+1]],0],{k,-kmax,kmax,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}];
  M1\[GothicR]=Table[If[-kmax<=k-kk<=kmax,If[k==kk,1,0]*((-I*kk*\[CapitalUpsilon]\[Theta]hat)^2+\[ScriptCapitalH]\[GothicR]func[p,0,0,Sqrt[Qhat]]),0],{k,-kmax,kmax,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}];
  M1\[CapitalUpsilon]\[Theta]=\[ScriptCapitalI]1\[Theta]four[[1;;-1;;2]];
  M1t=\[ScriptCapitalI]2four[[1;;-1;;2]];
  M1\[Phi]=\[ScriptCapitalI]3four[[1;;-1;;2]];
  M1=Transpose[Join[Transpose[-Im[M1\[Theta]]],Transpose[Re[M1\[GothicR]]],Re[{M1\[CapitalUpsilon]\[Theta],M1t,M1\[Phi]}]]];
  M2\[Theta]=Table[If[-kmax<=k-kk<=kmax,-kk^2*\[CapitalUpsilon]\[Theta]hat^2*\[ScriptCapitalQ]\[Theta]four[[k-kk+kmax+1]]-I*kk*\[CapitalUpsilon]\[Theta]hat*\[ScriptCapitalS]\[Theta]four[[k-kk+kmax+1]]+\[ScriptCapitalT]\[Theta]four[[k-kk+kmax+1]],0],{k,-kmax+1,kmax-1,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}];
  M2\[GothicR]=Table[0,{k,-kmax+1,kmax-1,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}];
  M2\[CapitalUpsilon]\[Theta]=\[ScriptCapitalU]1\[Theta]four[[2;;-2;;2]];
  M2t=\[ScriptCapitalU]2four[[2;;-2;;2]];
  M2\[Phi]=\[ScriptCapitalU]3four[[2;;-2;;2]];
  M2=Transpose[Join[Transpose[-Im@M2\[Theta]],Transpose[M2\[GothicR]],Re@{M2\[CapitalUpsilon]\[Theta],M2t,M2\[Phi]}]];
  M3\[Theta]=Table[If[-kmax<=k-kk<=kmax,-I*kk*\[CapitalUpsilon]\[Theta]hat*\[ScriptCapitalK]\[Theta]four[[k-kk+kmax+1]]+\[ScriptCapitalM]\[Theta]four[[k-kk+kmax+1]],0],{k,-kmax,kmax,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}];
  M3\[GothicR]=Table[0,{k,-kmax,kmax,2},{kk,Join[Range[-kmax,-2,2],Range[2,kmax,2]]}];
  M3\[CapitalUpsilon]\[Theta]=\[ScriptCapitalN]1\[Theta]four[[1;;-1;;2]];
  M3t=\[ScriptCapitalN]2four[[1;;-1;;2]];
  M3\[Phi]=\[ScriptCapitalN]3four[[1;;-1;;2]];
  M3=Transpose[Join[Transpose[-Im@M3\[Theta]],Transpose[M3\[GothicR]],Re@{M3\[CapitalUpsilon]\[Theta],M3t,M3\[Phi]}]];
  M=Join[M1,M2,M3];
  c=-Re[Join[Flatten[\[ScriptCapitalJ]four[[1;;-1;;2]]],Flatten[\[ScriptCapitalV]four[[2;;-2;;2]]],Flatten[\[ScriptCapitalP]four[[1;;-1;;2]]]]];(* source term, even, odd, even k-modes for r, \[Theta], norm equations *)
  v=LeastSquares[M,c];
  \[Delta]\[Chi]\[Theta]SfourIm=v[[1;;kmax]];
  \[Delta]\[GothicR]four=v[[kmax+1;;kmax+kmax]];
  \[CapitalUpsilon]\[Theta]S=Re[v[[-3]]];
  udt0S=Re[v[[-2]]];
  ud\[Phi]0S=Re[v[[-1]]];
  <|
    "\[ScriptCapitalR]tfourIm"->Im[\[ScriptCapitalR]tfour],
    "\[ScriptCapitalR]\[Phi]fourIm"->Im[\[ScriptCapitalR]\[Phi]four],
    "\[Delta]\[Chi]\[Theta]SfourIm"->\[Delta]\[Chi]\[Theta]SfourIm,
    "\[Delta]\[GothicR]four"->\[Delta]\[GothicR]four,
    "\[CapitalUpsilon]\[Theta]S"->\[CapitalUpsilon]\[Theta]S,
    "udt0S"->udt0S,
    "ud\[Phi]0S"->ud\[Phi]0S,
    "error"->Norm[M . v-c],(* error estimate from the residual *)
    "kmax"->kmax,
    "orbitGeo"->orbitGeo
  |>
]


KerrSpinOrbitCorrection2Spherical[correction_]:=Module[{a,p,e,x,\[ScriptCapitalI],prec,\[ScriptCapitalR]\[Phi]fourIm,\[ScriptCapitalR]tfourIm,kmax,orbit,Ehat,Lzhat,Qhat,Khat,\[CapitalUpsilon]\[Theta]hat,\[CapitalUpsilon]\[Phi]hat,\[CapitalGamma]hat,
    \[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]\[Theta],z1,z2,kz,steps\[Theta],\[Delta]udtS,\[Delta]ud\[Phi]S,\[Delta]\[Chi]\[Theta]SfourIm,\[Delta]\[Chi]\[Theta]S,d\[Delta]\[Chi]\[Theta]S,\[Delta]\[GothicR]four,\[Delta]\[GothicR],d\[Delta]\[GothicR],\[CapitalUpsilon]\[Theta]S,\[CapitalUpsilon]\[Phi]S,\[CapitalGamma]S,
    udt0S,ud\[Phi]0S,dzd\[Lambda],d\[Theta]d\[Lambda],zhat,\[Chi]\[Theta]hat,St,Sr,Sz,S\[Phi],\[CapitalDelta],\[CapitalSigma],rS,zS,UtS,U\[Phi]S,UtSlist,U\[Phi]Slist,UtSfour,U\[Phi]Sfour,Im\[CapitalDelta]tSfour,Im\[CapitalDelta]\[Phi]Sfour,\[CapitalDelta]tS,\[CapitalDelta]\[Phi]S,Uthat,U\[Phi]hat,
    Uthatlist,U\[Phi]hatlist,Uthatfour,U\[Phi]hatfour,\[CapitalSigma]gtt,\[CapitalSigma]g\[Phi]\[Phi],\[CapitalSigma]gt\[Phi],d\[CapitalSigma]gttdr,d\[CapitalSigma]gttdz,d\[CapitalSigma]g\[Phi]\[Phi]dr,d\[CapitalSigma]g\[Phi]\[Phi]dz,d\[CapitalSigma]gt\[Phi]dr,d\[CapitalSigma]gt\[Phi]dz,Kkz,h1,h2,h3,ES,LS,OrbitCorrection,
    ExpjkTable,udtS,ud\[Phi]S},
  h1[r_,z_]:=(r (-3 a^2 r^2 z^2+a^4 z^4+Khat (r^2-3 a^2 z^2)))/(Sqrt[Khat] (r^2+a^2 z^2)^3);(* functions for the four-velocity equations *)
  h2[r_,z_]:=1/(Sqrt[Khat] (r^2+a^2 z^2)^3) (-Ehat Lzhat r^6+a^4 Ehat Lzhat r^2 z^4+a^2 Ehat Lzhat r^4 (-2+z^2)-a^6 Ehat Lzhat z^4 (-2+z^2)+a^7 Ehat^2 z^4 (-1+z^2)+a r^3 (Lzhat^2 r+Khat (-1+z^2)+Khat r (-1+2 z^2)+r^3 (z^2-Ehat^2 (-1+z^2)))+a^3 (-Ehat^2 r^4 (-1+z^2)+r z^2 (2 r^3+2 Khat r z^2-3 Khat (-1+z^2)-3 r^2 (-1+z^2)))+a^5 z^4 (Khat-Lzhat^2+r ((-1+z^2)+r (2-z^2+Ehat^2 (-1+z^2)))));
  h3[r_,z_]:=-((2 a r z)/(Sqrt[Khat] (r^2+a^2 z^2)^2));
  \[Chi]\[Theta]hat[qz_]:=Floor[qz/(2Pi)]*2Pi+If[Mod[qz,2Pi]<Pi,ArcCos[JacobiSN[Kkz*2*(qz+Pi/2)/Pi,kz]],2Pi-ArcCos[JacobiSN[Kkz*2*(qz+Pi/2)/Pi,kz]]];
  \[Delta]\[Chi]\[Theta]S[wz_]:=2*Sum[\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*Sin[2k*wz],{k,1,kmax/2}];
  d\[Delta]\[Chi]\[Theta]S[wz_]:=2*Sum[2*k*\[CapitalUpsilon]\[Theta]hat*\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*Cos[2k*wz],{k,1,kmax/2}];
  \[Delta]\[GothicR][wz_]:=2*Sum[\[Delta]\[GothicR]four[[k+kmax/2]]*Cos[2k*wz],{k,1,kmax/2}];
  d\[Delta]\[GothicR][wz_]:=-2*Sum[(2k*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicR]four[[k+kmax/2]]*Sin[2k*wz],{k,1,kmax/2}];
  dzd\[Lambda][wz_]:=-Sign[Pi-Mod[wz,2Pi]]*Sqrt[Qhat-zhat[wz]^2*(a^2*(1-Ehat^2)*(1-zhat[wz]^2)+Lzhat^2+Qhat)];
  d\[Theta]d\[Lambda][wz_]:=-dzd\[Lambda][wz]/Sqrt[1-zhat[wz]^2];
  zhat[wz_]:=Cos[orbit["Trajectory"][[3]][wz]];
  St[wz_]:=-1/Sqrt[Khat]*(-a*(p*dzd\[Lambda][wz]+zhat[wz]*0)/\[CapitalSigma][p,zhat[wz]]);
  Sr[wz_]:=-1/Sqrt[Khat]*(a*zhat[wz]*((p^2+a^2)*Ehat-a*Lzhat)/\[CapitalDelta][p]);
  Sz[wz_]:=-1/Sqrt[Khat]*(a*p*Ehat-(p*Lzhat)/(1-zhat[wz]^2));
  S\[Phi][wz_]:=-1/Sqrt[Khat]*((a^2*zhat[wz]*(1-zhat[wz]^2)*0+p*(p^2+a^2)*dzd\[Lambda][wz])/\[CapitalSigma][p,zhat[wz]]);
  \[CapitalDelta][r_]:=r^2-2*r+a^2;
  \[CapitalSigma][r_,z_]:=r^2+a^2*z^2;
  \[CapitalSigma]gtt[r_,z_]:=-((r^2+a^2)^2/(a^2+(-2+r) r))+a^2 (1-z^2);
  \[CapitalSigma]g\[Phi]\[Phi][r_,z_]:=1/(1-z^2)-a^2/(a^2+r (-2+r));
  \[CapitalSigma]gt\[Phi][r_]:=a(1-(r^2+a^2)/(a^2-2 r+r^2));
  d\[CapitalSigma]gttdr[r_]:=2 ((r^2+a^2) ((r^2+a^2) (r-1)-2 (r^2-2r+a^2) r))/(r^2-2r+a^2)^2;
  d\[CapitalSigma]gttdz[z_]:=-2 a^2 z;
  d\[CapitalSigma]gt\[Phi]dr[r_]:=2a ((r^2+a^2) (r-1)-(r^2-2r+a^2) r)/(r^2-2r+a^2)^2;
  d\[CapitalSigma]g\[Phi]\[Phi]dr[r_]:=(2 a^2 (r-1))/(a^2+(-2+r) r)^2;
  d\[CapitalSigma]g\[Phi]\[Phi]dz[z_]:=(2 z)/(-1+z^2)^2;
  rS[wz_]:=\[Delta]\[GothicR][wz];
  zS[wz_]:=-Sin[\[ScriptCapitalI]]*\[Delta]\[Chi]\[Theta]S[wz]*Sin[\[Chi]\[Theta]hat[wz]];
  Uthat[wz_]:=\[CapitalSigma]gtt[p,zhat[wz]]*(-Ehat )+\[CapitalSigma]gt\[Phi][p]*(Lzhat);
  U\[Phi]hat[wz_]:=\[CapitalSigma]g\[Phi]\[Phi][p,zhat[wz]]*(Lzhat)+\[CapitalSigma]gt\[Phi][p]*(-Ehat);
  \[CapitalDelta]tS[wz_]:=Sum[2*Im\[CapitalDelta]tSfour[[k+kmax+1]]*Sin[k*wz],{k,1,kmax}]+Sum[2*Im\[CapitalDelta]tSfour[[k+kmax+1]]*Sin[k*wz],{k,-kmax,kmax}];
  \[CapitalDelta]\[Phi]S[wz_]:=Sum[2*Im\[CapitalDelta]\[Phi]Sfour[[k+kmax+1]]*Sin[k*wz],{k,1,kmax}]+Sum[2*Im\[CapitalDelta]\[Phi]Sfour[[k+kmax+1]]*Sin[k*wz],{k,-kmax,kmax}];
  OrbitCorrection[wz_]:=Module[{sinkwz,coskwz,\[Delta]\[Chi]\[Theta]Sp,\[Chi]\[Theta]hatp,k},
    sinkwz={Sin[wz]};
    coskwz={Cos[wz]};
    For[k=2,k<=kmax,k++,
      AppendTo[sinkwz,sinkwz[[k-1]]*coskwz[[1]]+coskwz[[k-1]]*sinkwz[[1]]];
      AppendTo[coskwz,coskwz[[k-1]]*coskwz[[1]]-sinkwz[[k-1]]*sinkwz[[1]]];
    ];
    Clear[k];
    \[Delta]\[Chi]\[Theta]Sp=(2*Sum[\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*sinkwz[[2k]],{k,1,kmax/2}]);
    \[Chi]\[Theta]hatp=\[Chi]\[Theta]hat[wz];
    <|
      "\[Delta]\[CapitalDelta]t"->2*(Sum[Im\[CapitalDelta]tSfour[[k+kmax+1]]*sinkwz[[k]],{k,1,kmax}]),
      "\[Delta]r"->2*(Sum[\[Delta]\[GothicR]four[[k+kmax/2]]*coskwz[[2k]],{k,1,kmax/2}]),
      "\[Delta]z"->-Sin[\[ScriptCapitalI]]*\[Delta]\[Chi]\[Theta]Sp*Sin[\[Chi]\[Theta]hatp],
      "\[Delta]\[CapitalDelta]\[Phi]"->2*(Sum[Im\[CapitalDelta]\[Phi]Sfour[[k+kmax+1]]*sinkwz[[k]],{k,1,kmax}]),
      "\[Delta]Ur"->(-2*(Sum[(2k*\[CapitalUpsilon]\[Theta]hat)*\[Delta]\[GothicR]four[[k+kmax/2]]*sinkwz[[2k]],{k,1,kmax/2}])),
      "\[Delta]Uz"->(-Sin[\[ScriptCapitalI]]*(\[Chi]\[Theta]hat'[wz]*(Cos[\[Chi]\[Theta]hatp]*\[Delta]\[Chi]\[Theta]Sp*\[CapitalUpsilon]\[Theta]hat+Sin[\[Chi]\[Theta]hatp]*\[CapitalUpsilon]\[Theta]S)+Sin[\[Chi]\[Theta]hatp]*(2*Sum[2*k*\[CapitalUpsilon]\[Theta]hat*\[Delta]\[Chi]\[Theta]SfourIm[[k+kmax/2]]*coskwz[[2k]],{k,1,kmax/2}])))
    |>
  ];
  UtS[wz_]:=((d\[CapitalSigma]gttdr[p]*rS[wz]+d\[CapitalSigma]gttdz[zhat[wz]]*zS[wz])*(-Ehat)+\[CapitalSigma]gtt[p,zhat[wz]]*(udtS[wz])+d\[CapitalSigma]gt\[Phi]dr[p]*rS[wz]*(Lzhat)+\[CapitalSigma]gt\[Phi][p]*(ud\[Phi]S[wz]));
  U\[Phi]S[wz_]:=((d\[CapitalSigma]g\[Phi]\[Phi]dr[p]*rS[wz]+d\[CapitalSigma]g\[Phi]\[Phi]dz[zhat[wz]]*zS[wz])*(Lzhat)+\[CapitalSigma]g\[Phi]\[Phi][p,zhat[wz]]*(ud\[Phi]S[wz])+d\[CapitalSigma]gt\[Phi]dr[p]*rS[wz]*(-Ehat)+\[CapitalSigma]gt\[Phi][p]*(udtS[wz]));
(****************************************)
  orbit=correction["orbitGeo"];
  a=orbit["a"];
  p=orbit["p"];
  e=orbit["e"];
  x=orbit["Inclination"];
  \[ScriptCapitalI]=ArcCos[x];
  prec=Precision[{a,p,x}];
  {Ehat,Lzhat,Qhat}=Values[orbit["ConstantsOfMotion"]];
  Khat=Qhat+(a*Ehat-Lzhat)^2;
  \[CapitalUpsilon]\[Theta]hat=orbit["PolarFrequency"];
  \[CapitalUpsilon]\[Phi]hat=orbit["AzimuthalFrequency"];
  \[CapitalGamma]hat=orbit["Frequencies"]["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"];
  \[CapitalDelta]t\[Theta]=orbit["TrajectoryDeltas"]["\[CapitalDelta]t\[Theta]"];
  \[CapitalDelta]\[Phi]\[Theta]=orbit["TrajectoryDeltas"]["\[CapitalDelta]\[Phi]\[Theta]"];
  z1=Sqrt[1-x^2];
  z2=Sqrt[a^2*(1-Ehat^2)+Lzhat^2/(1-z1^2)];
  kz=a^2*(1-Ehat^2)*z1^2/z2^2;
  Kkz=EllipticK[kz];
  kmax=correction["kmax"];
  steps\[Theta]=4*kmax;
  \[ScriptCapitalR]\[Phi]fourIm=correction["\[ScriptCapitalR]\[Phi]fourIm"];
  \[ScriptCapitalR]tfourIm=correction["\[ScriptCapitalR]tfourIm"];
  \[Delta]\[Chi]\[Theta]SfourIm=correction["\[Delta]\[Chi]\[Theta]SfourIm"];
  \[Delta]\[GothicR]four=correction["\[Delta]\[GothicR]four"];
  \[CapitalUpsilon]\[Theta]S=correction["\[CapitalUpsilon]\[Theta]S"];
  udt0S=correction["udt0S"];
  ud\[Phi]0S=correction["ud\[Phi]0S"];
  Print["Calculating ES,LS"];
  ES = -udt0S + Sum[h1[p,zhat[w\[Theta]]], {w\[Theta],2Pi/steps\[Theta]/2,2Pi,2Pi/steps\[Theta]}]/steps\[Theta];
  LS = ud\[Phi]0S - Sum[h2[p,zhat[w\[Theta]]], {w\[Theta],2Pi/steps\[Theta]/2,2Pi,2Pi/steps\[Theta]}]/steps\[Theta];
  udtS[wz_]:=-(ES)+h1[p,zhat[wz]];
  ud\[Phi]S[wz_]:=(LS)+(h2[p,zhat[wz]]);
  Print["Calculating UtSlist, U\[Phi]Slist"];
  UtSlist=Table[UtS[wz],{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  U\[Phi]Slist=Table[U\[Phi]S[wz],{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  Print["Calculating Uthatlist, U\[Phi]hatlist"];
  Uthatlist=Table[Uthat[wz],{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  U\[Phi]hatlist=Table[U\[Phi]hat[wz],{wz,2Pi/(steps\[Theta])/2,2Pi,2Pi/(steps\[Theta])}];
  Print["Calculating UtSfour, U\[Phi]Sfour"];
  ExpjkTable=Table[N[Exp[2Pi*I*k*(j-1/2)/steps\[Theta]],Precision[{a,p,e,x}]],{j,1,steps\[Theta]},{k,-kmax,kmax}];
  UtSfour=UtSlist . ExpjkTable/(steps\[Theta]);
  U\[Phi]Sfour=U\[Phi]Slist . ExpjkTable/(steps\[Theta]);
  Print["Calculating Uthatfour, U\[Phi]hatfour"];
  Uthatfour=Uthatlist . ExpjkTable/(steps\[Theta]);
  U\[Phi]hatfour=U\[Phi]hatlist . ExpjkTable/(steps\[Theta]);
  Print["Calculating Im\[CapitalDelta]tSfour, Im\[CapitalDelta]\[Phi]Sfour"];
  Im\[CapitalDelta]tSfour=Table[If[k==0,0,Re@UtSfour[[k+kmax+1]]/(k*\[CapitalUpsilon]\[Theta]hat)-Re@Uthatfour[[k+kmax+1]]/(k*\[CapitalUpsilon]\[Theta]hat)^2*(k*\[CapitalUpsilon]\[Theta]S)],{k,-kmax,kmax}];
  Im\[CapitalDelta]\[Phi]Sfour=Table[If[k==0,0,Re@U\[Phi]Sfour[[k+kmax+1]]/(k*\[CapitalUpsilon]\[Theta]hat)-Re@U\[Phi]hatfour[[k+kmax+1]]/(k*\[CapitalUpsilon]\[Theta]hat)^2*(k*\[CapitalUpsilon]\[Theta]S)],{k,-kmax,kmax}];
  \[CapitalGamma]S=Re[UtSfour[[kmax+1]]];
  \[CapitalUpsilon]\[Phi]S=Re[U\[Phi]Sfour[[kmax+1]]];
  <|
    "a"->a,
    "p"->p,
    "e"->e,
    "Inclination"->x,
    "En0"->Ehat,
    "Lz0"->Lzhat,
    "K0"->Khat,
    "TrajectoryCorrection"->{Function[{wz},\[CapitalDelta]tS[wz]],Function[{wz},rS[wz]],Function[{wz},zS[wz]],Function[{wz},\[CapitalDelta]\[Phi]S[wz]]},
    "TrajectoryGeo"->orbit["Trajectory"],
    "MinoFrequenciesGeo"->{\[CapitalGamma]hat,\[CapitalUpsilon]\[Theta]hat,\[CapitalUpsilon]\[Phi]hat},
    "BLFrequenciesGeo"->{\[CapitalUpsilon]\[Theta]hat/\[CapitalGamma]hat,\[CapitalUpsilon]\[Phi]hat/\[CapitalGamma]hat},
    "FourVelocityCorrection"->{
      Function[{wz},(udtS[wz])],
      Function[{wz},d\[Delta]\[GothicR][wz]],
      Function[{wz},(-Sin[\[ScriptCapitalI]]*(\[Chi]\[Theta]hat'[wz]*(Cos[\[Chi]\[Theta]hat[wz]]*\[Delta]\[Chi]\[Theta]S[wz]*\[CapitalUpsilon]\[Theta]hat+Sin[\[Chi]\[Theta]hat[wz]]*\[CapitalUpsilon]\[Theta]S)+Sin[\[Chi]\[Theta]hat[wz]]*d\[Delta]\[Chi]\[Theta]S[wz]))],
      Function[{wz},(ud\[Phi]S[wz])]},
    "Spin4vector"->{Function[{wz},St[wz]],Function[{wz},Sr[wz]],Function[{wz},Sz[wz]],Function[{wz},S\[Phi][wz]]},
    "TrajectoryDeltas"->orbit["TrajectoryDeltas"],
    "MinoFrequenciesCorrection"->{\[CapitalGamma]S,\[CapitalUpsilon]\[Theta]S,\[CapitalUpsilon]\[Phi]S},
    "BLFrequenciesCorrection"->{\[CapitalUpsilon]\[Theta]S/\[CapitalGamma]hat-\[CapitalUpsilon]\[Theta]hat/\[CapitalGamma]hat^2*\[CapitalGamma]S,\[CapitalUpsilon]\[Phi]S/\[CapitalGamma]hat-\[CapitalUpsilon]\[Phi]hat/\[CapitalGamma]hat^2*\[CapitalGamma]S},
    "\[Delta]En"->ES,
    "\[Delta]Lz"->LS,
    "OrbitCorrection"->Function[{wz},OrbitCorrection[wz]],
    "error"->correction["error"]
  |>
]


KerrSpinOrbitCorrectionSpherical[a_, p_, e_, x_, kmax_?EvenQ]:=Module[{orbit,correction},
  orbit = KerrGeodesics`KerrGeoOrbit`KerrGeoOrbit[a, p, e, x, "Parametrization"->"Phases", "Method"->"Analytic"];
  correction = CorrectionSpherical[orbit, kmax];
  KerrSpinOrbitCorrection2Spherical[correction]
]


KerrSpinOrbitSpherical[orbitCorrection_,spar_]:=Module[{\[CapitalGamma],\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]},
  {\[CapitalGamma],\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]}=orbitCorrection["MinoFrequenciesGeo"]+spar*orbitCorrection["MinoFrequenciesCorrection"];
  \[CapitalOmega]\[Theta]=\[CapitalUpsilon]\[Theta]/\[CapitalGamma];
  \[CapitalOmega]\[Phi]=\[CapitalUpsilon]\[Phi]/\[CapitalGamma];
  <|
    "a"->orbitCorrection["a"],
    "p"->orbitCorrection["p"],
    "e"->orbitCorrection["e"],
    "Inclination"->orbitCorrection["Inclination"],
    "s"->spar,
    "En0"->orbitCorrection["En0"],
    "Lz0"->orbitCorrection["Lz0"],
    "K0"->orbitCorrection["K0"],
    "BLFrequencies"->{\[CapitalOmega]\[Theta],\[CapitalOmega]\[Phi]},
    "TrajectoryGeo"->orbitCorrection["TrajectoryGeo"],
    "TrajectoryDeltas"->orbitCorrection["TrajectoryDeltas"],
    "MinoFrequencies"->{\[CapitalGamma],\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]},
    "MinoFrequenciesGeo"->orbitCorrection["MinoFrequenciesGeo"],
    "MinoFrequenciesCorrection"->orbitCorrection["MinoFrequenciesCorrection"],
    "\[Delta]En"->orbitCorrection["\[Delta]En"],
    "LS"->orbitCorrection["\[Delta]Lz"],
    "OrbitCorrection"->orbitCorrection["OrbitCorrection"]
  |>
]


(* ::Subsection::Closed:: *)
(*Analytical equatorial orbit*)


KerrSpinOrbitEquatorial[orbit_]:=Module[{M=1,a,En,L,Q,r1,r2,r3,r4,kr,rp,rm,hr,hp,hm,traj},
  If[Abs[orbit["Inclination"]]!=1&,Print["Cannot use offequatorial orbit as input"];Return[$Failed]];
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


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];

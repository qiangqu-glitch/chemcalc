import { useState, useMemo, useCallback, useRef, useEffect } from "react";
import { XAxis, YAxis, Tooltip, ResponsiveContainer, LineChart, Line, CartesianGrid, ReferenceLine } from "recharts";

var R_GAS = 8.314462618;

// Component data loaded at runtime from API Worker
// Data endpoint: https://api.chemcalc.cn/api/data
var API_URL = "https://api.chemcalc.cn/api/data";
var COMPS = [];
var KIJ = {};
var LIQDEN = {};
var _dataLoaded = false;
var _dataLoadPromise = null;

function loadChemData() {
  if (_dataLoaded) return Promise.resolve();
  if (_dataLoadPromise) return _dataLoadPromise;
  _dataLoadPromise = fetch(API_URL)
    .then(function(r){ return r.json(); })
    .then(function(d){
      COMPS = d.COMPS || [];
      KIJ = d.KIJ || {};
      LIQDEN = d.LIQDEN || {};
      _dataLoaded = true;
    })
    .catch(function(e){
      console.error("ChemCalc: failed to load component data", e);
    });
  return _dataLoadPromise;
}


// Liquid density (kg/m3) from DIPPR 105 or Rackett fallback
function liquidDensity(T, comp) {
  if (T >= comp.Tc) return null;
  var dd = LIQDEN[comp.id];
  if (dd) {
    // DIPPR 105: rho = A / B^(1+(1-T/C)^D) [kmol/m3]; kmol/m3 * g/mol = kg/m3 (numerically exact)
    var rhoMol = dd.A / Math.pow(dd.B, 1 + Math.pow(1-T/dd.C, dd.D));
    return rhoMol * comp.MW; // kg/m3
  }
  // Rackett fallback: V = R*Tc/Pc * ZRA^(1+(1-Tr)^(2/7))
  var Tr = T / comp.Tc;
  var Zra = 0.29056 - 0.08775 * comp.w;
  var Vsat = R_GAS * comp.Tc / comp.Pc * Math.pow(Zra, 1 + Math.pow(1-Tr, 2/7));
  return comp.MW / 1000 / Vsat;
}

function getKij(a, b) { return a === b ? 0 : (KIJ[a+"-"+b] || KIJ[b+"-"+a] || 0); }

var TU = [
  { id:"K", lb:"K", toK:function(v){return v}, fromK:function(v){return v} },
  { id:"C", lb:"\u00b0C", toK:function(v){return v+273.15}, fromK:function(v){return v-273.15} },
  { id:"F", lb:"\u00b0F", toK:function(v){return (v-32)*5/9+273.15}, fromK:function(v){return (v-273.15)*9/5+32} },
];
var PU = [
  { id:"MPa", lb:"MPa(a)", toP:function(v){return v*1e6}, fromP:function(v){return v/1e6} },
  { id:"MPag", lb:"MPa(g)", toP:function(v){return (v+0.101325)*1e6}, fromP:function(v){return v/1e6-0.101325} },
  { id:"kPa", lb:"kPa(a)", toP:function(v){return v*1e3}, fromP:function(v){return v/1e3} },
  { id:"kPag", lb:"kPa(g)", toP:function(v){return (v+101.325)*1e3}, fromP:function(v){return v/1e3-101.325} },
  { id:"bar", lb:"bar(a)", toP:function(v){return v*1e5}, fromP:function(v){return v/1e5} },
  { id:"barg", lb:"bar(g)", toP:function(v){return (v+1.01325)*1e5}, fromP:function(v){return v/1e5-1.01325} },
  { id:"atm", lb:"atm(a)", toP:function(v){return v*101325}, fromP:function(v){return v/101325} },
  { id:"psi", lb:"psia", toP:function(v){return v*6894.757}, fromP:function(v){return v/6894.757} },
  { id:"psig", lb:"psig", toP:function(v){return (v+14.696)*6894.757}, fromP:function(v){return v/6894.757-14.696} },
  { id:"mmHg", lb:"mmHg(a)", toP:function(v){return v*133.322}, fromP:function(v){return v/133.322} },
];

function prAB(T, Tc, Pc, w) {
  var Tr = T / Tc, m = 0.37464 + 1.54226*w - 0.26992*w*w;
  var s = 1 + m*(1 - Math.sqrt(Tr)), alpha = s*s;
  return { a: 0.45724 * R_GAS*R_GAS * Tc*Tc / Pc * alpha, b: 0.07780 * R_GAS * Tc / Pc };
}

function srkAB(T, Tc, Pc, w) {
  var Tr = T / Tc, m = 0.480 + 1.574*w - 0.176*w*w;
  var s = 1 + m*(1 - Math.sqrt(Tr)), alpha = s*s;
  return { a: 0.42748 * R_GAS*R_GAS * Tc*Tc / Pc * alpha, b: 0.08664 * R_GAS * Tc / Pc };
}

function solveCubic(c2, c1, c0) {
  var p = (3*c1 - c2*c2) / 3, q = (2*c2*c2*c2 - 9*c2*c1 + 27*c0) / 27;
  var disc = q*q/4 + p*p*p/27, roots = [];
  if (disc > 1e-14) {
    var sd = Math.sqrt(disc);
    roots.push(Math.cbrt(-q/2 + sd) + Math.cbrt(-q/2 - sd) - c2/3);
  } else {
    var r = Math.sqrt(Math.max(0, -p*p*p/27));
    var theta = Math.acos(Math.max(-1, Math.min(1, -q/(2*Math.max(r,1e-30)))));
    var mm = 2 * Math.cbrt(Math.max(r,0));
    roots.push(mm*Math.cos(theta/3) - c2/3);
    roots.push(mm*Math.cos((theta+2*Math.PI)/3) - c2/3);
    roots.push(mm*Math.cos((theta+4*Math.PI)/3) - c2/3);
  }
  return roots.filter(function(z){return z > 0}).sort(function(a,b){return a-b});
}

function prRoots(A, B) { return solveCubic(-(1-B), A-3*B*B-2*B, -(A*B-B*B-B*B*B)); }
function srkRoots(A, B) { return solveCubic(-1, A-B-B*B, -A*B); }

function prLnPhi(Z, A, B) {
  var s2 = 1.4142135623730951, n = Z + (1+s2)*B, d = Z + (1-s2)*B;
  if (n <= 0 || d <= 0 || Z-B <= 0) return NaN;
  return (Z-1) - Math.log(Z-B) - A/(2*s2*B)*Math.log(n/d);
}
function srkLnPhi(Z, A, B) {
  if (Z-B <= 0 || Z+B <= 0 || Z <= 0) return NaN;
  return (Z-1) - Math.log(Z-B) - A/B*Math.log((Z+B)/Z);
}

function calcZ(T, P, comp, phase, method) {
  var ab = method === "SRK" ? srkAB(T,comp.Tc,comp.Pc,comp.w) : prAB(T,comp.Tc,comp.Pc,comp.w);
  var A = ab.a*P/(R_GAS*R_GAS*T*T), B = ab.b*P/(R_GAS*T);
  var roots = method === "SRK" ? srkRoots(A,B) : prRoots(A,B);
  if (roots.length === 0) return null;
  return phase === "liquid" ? roots[0] : roots[roots.length-1];
}

function calcPsat(T, comp, method) {
  if (T >= comp.Tc*0.999 || T < 0.3*comp.Tc) return null;
  var Tr = T/comp.Tc;
  var f0 = 5.92714 - 6.09648/Tr - 1.28862*Math.log(Tr) + 0.169347*Math.pow(Tr,6);
  var f1 = 15.2518 - 15.6875/Tr - 13.4721*Math.log(Tr) + 0.43577*Math.pow(Tr,6);
  var P = comp.Pc * Math.exp(f0 + comp.w*f1);
  if (isNaN(P) || P <= 0 || P >= comp.Pc) P = comp.Pc*Math.exp(5.373*(1+comp.w)*(1-comp.Tc/T));
  P = Math.max(P, 100); P = Math.min(P, comp.Pc*0.95);
  var getAB = method === "SRK" ? srkAB : prAB;
  var getRoots = method === "SRK" ? srkRoots : prRoots;
  var getLnPhi = method === "SRK" ? srkLnPhi : prLnPhi;
  var ab = getAB(T, comp.Tc, comp.Pc, comp.w);
  for (var it = 0; it < 200; it++) {
    var A = ab.a*P/(R_GAS*R_GAS*T*T), B = ab.b*P/(R_GAS*T);
    var roots = getRoots(A, B);
    if (roots.length < 2) { P = P > comp.Pc*0.5 ? P*0.7 : P*1.5; P = Math.max(100, Math.min(P, comp.Pc*0.95)); continue; }
    var Zl = roots[0], Zv = roots[roots.length-1];
    if (Math.abs(Zv-Zl) < 1e-8) return P;
    var lpL = getLnPhi(Zl, A, B), lpV = getLnPhi(Zv, A, B);
    if (isNaN(lpL) || isNaN(lpV)) { P *= 0.9; continue; }
    var diff = lpL - lpV;
    if (Math.abs(diff) < 1e-9) return P;
    var P2 = Math.exp(Math.log(P) + diff/Math.max(Zv-Zl, 0.01));
    if (P2/P > 3) P2 = P*3; if (P2/P < 0.33) P2 = P*0.33;
    P = Math.max(10, Math.min(P2, comp.Pc*0.999));
  }
  return P;
}

function cpIg(T, c) { return R_GAS*(c.cp[0]+c.cp[1]*T+c.cp[2]*T*T+c.cp[3]*T*T*T+c.cp[4]*Math.pow(T,4)); }
function hIg(T, c) { var T0=298.15, p=c.cp; return R_GAS*(p[0]*(T-T0)+p[1]/2*(T*T-T0*T0)+p[2]/3*(T*T*T-T0*T0*T0)+p[3]/4*(Math.pow(T,4)-Math.pow(T0,4))+p[4]/5*(Math.pow(T,5)-Math.pow(T0,5))); }
function sIg(T, P, c) { var T0=298.15, P0=101325, p=c.cp; return R_GAS*(p[0]*Math.log(T/T0)+p[1]*(T-T0)+p[2]/2*(T*T-T0*T0)+p[3]/3*(T*T*T-T0*T0*T0)+p[4]/4*(Math.pow(T,4)-Math.pow(T0,4))) - R_GAS*Math.log(P/P0); }

function prDep(T, P, Z, comp) {
  var ab = prAB(T, comp.Tc, comp.Pc, comp.w);
  var A = ab.a*P/(R_GAS*R_GAS*T*T), B = ab.b*P/(R_GAS*T), s2 = 1.4142135623730951;
  var Tr = T/comp.Tc, m = 0.37464+1.54226*comp.w-0.26992*comp.w*comp.w;
  var sqAlpha = Math.sqrt(ab.a * comp.Pc / (0.45724*R_GAS*R_GAS*comp.Tc*comp.Tc));
  var dadT = -0.45724*R_GAS*R_GAS*comp.Tc*comp.Tc/comp.Pc*m*sqAlpha/(Math.sqrt(Tr)*comp.Tc);
  var n = Z+(1+s2)*B, d = Z+(1-s2)*B;
  if (n <= 0 || d <= 0) return { Hd:0, Sd:0 };
  var lnT = Math.log(n/d);
  return { Hd: R_GAS*T*(Z-1) + (T*dadT-ab.a)/(2*s2*ab.b)*lnT, Sd: R_GAS*Math.log(Math.max(Z-B,1e-30)) + dadT/(2*s2*ab.b)*lnT };
}

function calcPure(T, P, comp, method) {
  var res = { T:T, P:P, comp:comp, method:method };
  var Tr = T/comp.Tc, Pr = P/comp.Pc;
  if (Tr >= 1 && Pr >= 1) { res.ph = "supercritical"; res.phCN = "\u8d85\u4e34\u754c"; }
  else if (Tr >= 1) { res.ph = "vapor"; res.phCN = "\u6c14\u76f8"; }
  else {
    var Ps = calcPsat(T, comp, method); res.Psat = Ps;
    if (Ps) {
      if (Math.abs(P-Ps)/Ps < 0.005) { res.ph = "two-phase"; res.phCN = "\u4e24\u76f8\u533a"; }
      else if (P < Ps) { res.ph = "vapor"; res.phCN = "\u6c14\u76f8"; }
      else { res.ph = "liquid"; res.phCN = "\u6db2\u76f8"; }
    } else { res.ph = P < comp.Pc ? "vapor" : "liquid"; res.phCN = res.ph === "vapor" ? "\u6c14\u76f8" : "\u6db2\u76f8"; }
  }
  var pForZ = (res.ph === "liquid" || res.ph === "two-phase") ? "liquid" : "vapor";
  var Z = calcZ(T, P, comp, pForZ, method);
  if (!Z || Z <= 0) { res.err = "Cannot solve EOS"; return res; }
  res.Z = Z; res.Vm = Z*R_GAS*T/P; res.rho = comp.MW/1000/res.Vm;
  // For liquid phase: use DIPPR correlation for accurate density (cubic EOS gives 5-15% error for liquids)
  if (res.ph === "liquid" || res.ph === "two-phase") {
    var liqRho = liquidDensity(T, comp);
    if (liqRho && liqRho > 0) { res.rho = liqRho; res.Vm = comp.MW/1000/liqRho; res.liqDenSrc = "DIPPR"; }
  }
  res.CpIg = cpIg(T, comp); res.Hig = hIg(T, comp); res.Sig = sIg(T, P, comp);
  var dep = prDep(T, P, Z, comp);
  res.H = res.Hig + dep.Hd; res.S = res.Sig + dep.Sd;
  res.Cp = res.CpIg; res.Cv = res.Cp - R_GAS; res.gamma = res.Cp / Math.max(res.Cv, 0.01);
  if (res.Psat && T < comp.Tc) {
    var Zv = calcZ(T, res.Psat, comp, "vapor", method);
    var Zl = calcZ(T, res.Psat, comp, "liquid", method);
    if (Zv && Zl) { res.rhoV = comp.MW/1000/(Zv*R_GAS*T/res.Psat); }
    var satLiqRho = liquidDensity(T, comp);
    if (satLiqRho) { res.rhoL = satLiqRho; } else if (Zl) { res.rhoL = comp.MW/1000/(Zl*R_GAS*T/res.Psat); }
  }
  return res;
}

function steamPsat(T) {
  if (T < 273.15 || T > 647.096) return null;
  var n = [0.11670521452767E+04,-0.72421316703206E+06,-0.17073846940092E+02,0.12020824702470E+05,-0.32325550322333E+07,0.14915108613530E+02,-0.48232657361591E+04,0.40511340542057E+06,-0.23855557567849E+00,0.65017534844798E+03];
  var th = T + n[8]/(T-n[9]);
  var AA = th*th+n[0]*th+n[1], BB = n[2]*th*th+n[3]*th+n[4], CC = n[5]*th*th+n[6]*th+n[7];
  return Math.pow(2*CC/(-BB+Math.sqrt(BB*BB-4*AA*CC)),4)*1e6;
}

// IAPWS-IF97 Region 4 backward: Tsat from P
function iapws97_Tsat(P) {
  var Pmpa = P / 1e6;
  if (Pmpa < 611.213e-6 || Pmpa > 22.064) return null;
  var n = [0.11670521452767E+04,-0.72421316703206E+06,-0.17073846940092E+02,0.12020824702470E+05,-0.32325550322333E+07,0.14915108613530E+02,-0.48232657361591E+04,0.40511340542057E+06,-0.23855557567849E+00,0.65017534844798E+03];
  var beta = Math.pow(Pmpa, 0.25);
  var E = beta*beta + n[2]*beta + n[5];
  var F = n[0]*beta*beta + n[3]*beta + n[6];
  var G = n[1]*beta*beta + n[4]*beta + n[7];
  var D = 2*G / (-F - Math.sqrt(F*F - 4*E*G));
  return (n[9]+D-Math.sqrt((n[9]+D)*(n[9]+D)-4*(n[8]+n[9]*D)))/2;
}

function steamR1(T, P) {
  var Ps=16.53e6, Ts=1386, pi=P/Ps, tau=Ts/T, Rw=461.526;
  var I=[0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,5,8,8,21,23,29,30,31,32];
  var J=[-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3,17,-4,0,6,-5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41];
  var nn=[0.14632971213167E+00,-0.84548187169114E+00,-0.37563603672040E+01,0.33855169168385E+01,-0.95791963387872E+00,0.15772038513228E+00,-0.16616417199501E-01,0.81214629983568E-03,0.28319080123804E-03,-0.60706301565874E-03,-0.18990068218419E-01,-0.32529748770505E-01,-0.21841717175414E-01,-0.52838357969930E-04,-0.47184321073267E-03,-0.30001780793026E-03,0.47661393906987E-04,-0.44141845330846E-05,-0.72694996297594E-15,-0.31679644845054E-04,-0.28270797985312E-05,-0.85205128120103E-09,-0.22425281908000E-05,-0.65171222895601E-06,-0.14341729937924E-12,-0.40516996860117E-06,-0.12734301741682E-08,-0.17424871230634E-09,-0.68762131295531E-18,0.14478307828521E-19,0.26335781662795E-22,-0.11947622640071E-22,0.18228094581404E-23,-0.93537087292458E-25];
  var gm=0,gp=0,gt=0,gtt=0,gpp=0,gpt=0;
  for (var i=0;i<34;i++) {
    var pt=Math.pow(7.1-pi,I[i]), tt=Math.pow(tau-1.222,J[i]);
    gm += nn[i]*pt*tt;
    if (I[i] > 0) gp += -nn[i]*I[i]*Math.pow(7.1-pi,I[i]-1)*tt;
    if (I[i] > 1) gpp += nn[i]*I[i]*(I[i]-1)*Math.pow(7.1-pi,I[i]-2)*tt;
    gt += nn[i]*pt*J[i]*Math.pow(tau-1.222,J[i]-1);
    gtt += nn[i]*pt*J[i]*(J[i]-1)*Math.pow(tau-1.222,J[i]-2);
    if (I[i] > 0) gpt += -nn[i]*I[i]*Math.pow(7.1-pi,I[i]-1)*J[i]*Math.pow(tau-1.222,J[i]-1);
  }
  var v=Rw*T/P*pi*gp, h=Rw*T*tau*gt, s=Rw*(tau*gt-gm), cp=-Rw*tau*tau*gtt;
  var cvN=Math.pow(gp-tau*gpt,2)/gpp; var cv=cp+Rw*cvN;
  var wN=gp*gp/((gp-tau*gpt)*(gp-tau*gpt)/(tau*tau*gtt)-gpp);
  var w=Math.sqrt(Math.abs(Rw*T*wN));
  var u=h-P*v;
  return { v:v, h:h/1000, s:s/1000, cp:cp/1000, cv:cv/1000, rho:1/v, u:u/1000, w:w, ph:"liquid", phCN:"\u6db2\u6001\u6c34", rg:1 };
}

function steamR2(T, P) {
  var Ps=1e6, Ts=540, pi=P/Ps, tau=Ts/T, Rw=461.526;
  var J0=[0,1,-5,-4,-3,-2,-1,2,3];
  var n0=[-0.96927686500217E+01,0.10086655024209E+02,-0.56087911283020E-02,0.71452738081455E-01,-0.40710498223928E+00,0.14240819171444E+01,-0.43839511319450E+01,-0.28408632460772E+00,0.21268463753307E-01];
  var Ir=[1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,7,7,7,8,8,9,10,10,10,16,16,18,20,20,20,21,22,23,24,24,24];
  var Jr=[0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7,3,16,35,0,11,25,8,36,13,4,10,14,29,50,57,20,35,48,21,53,39,26,40,58];
  var nr=[-0.17731742473213E-02,-0.17834862292358E-01,-0.45996013696365E-01,-0.57581259083432E-01,-0.50325278727930E-01,-0.33032641670203E-04,-0.18948987516315E-03,-0.39392777243355E-02,-0.43797295650573E-01,-0.26674547914087E-04,0.20481737692309E-07,0.43870667284435E-06,-0.32277677238570E-04,-0.15033924542148E-02,-0.40668253562649E-04,-0.78847309559367E-09,0.12790717852285E-07,0.48225372718507E-06,0.22922076337661E-05,-0.16714766451061E-10,-0.21171472321355E-02,-0.23895741934104E+02,-0.59059564324270E-17,-0.12621808899101E-05,-0.38946842435739E-01,0.11256211360459E-10,-0.82311340897998E+01,0.19809712802088E-07,0.10406965210174E-18,-0.10234747095929E-12,-0.10018179379511E-08,-0.80882908646985E-10,0.10693031879409E+00,-0.33662250574171E+00,0.89185845355421E-24,0.30629316876232E-12,-0.42002467698208E-05,-0.59056029685639E-25,0.37826947613457E-05,-0.12768608934681E-14,0.73087610595061E-28,0.55414715350778E-16,-0.94369707241210E-06];
  var g0=Math.log(pi), g0t=0, g0tt=0, g0p=1/pi;
  var gr=0, grp=0, grt=0, grtt=0, grpp=0, grpt=0;
  for (var i=0;i<9;i++) { g0 += n0[i]*Math.pow(tau,J0[i]); g0t += n0[i]*J0[i]*Math.pow(tau,J0[i]-1); g0tt += n0[i]*J0[i]*(J0[i]-1)*Math.pow(tau,J0[i]-2); }
  for (var j=0;j<43;j++) {
    var pp=Math.pow(pi,Ir[j]), tt=Math.pow(tau-0.5,Jr[j]);
    gr += nr[j]*pp*tt;
    grp += nr[j]*Ir[j]*Math.pow(pi,Ir[j]-1)*tt;
    grpp += nr[j]*Ir[j]*(Ir[j]-1)*Math.pow(pi,Ir[j]-2)*tt;
    grt += nr[j]*pp*Jr[j]*Math.pow(tau-0.5,Jr[j]-1);
    grtt += nr[j]*pp*Jr[j]*(Jr[j]-1)*Math.pow(tau-0.5,Jr[j]-2);
    grpt += nr[j]*Ir[j]*Math.pow(pi,Ir[j]-1)*Jr[j]*Math.pow(tau-0.5,Jr[j]-1);
  }
  var v=Rw*T/P*pi*(g0p+grp), h=Rw*T*tau*(g0t+grt), s=Rw*(tau*(g0t+grt)-(g0+gr)), cp=-Rw*tau*tau*(g0tt+grtt);
  var cvN=Math.pow(1+pi*grp-tau*pi*grpt,2), cvD=1-pi*pi*grpp;
  var cv=-Rw*tau*tau*(g0tt+grtt)-Rw*cvN/cvD;
  var wN=Rw*T*(1+2*pi*grp+pi*pi*grp*grp);
  var wD=cvD+cvN/(tau*tau*(g0tt+grtt));
  var w=Math.sqrt(Math.abs(wN/wD));
  var u=h-P*v;
  return { v:v, h:h/1000, s:s/1000, cp:cp/1000, cv:cv/1000, rho:1/v, u:u/1000, w:w, ph:"vapor", phCN:"\u8fc7\u70ed\u84b8\u6c7d", rg:2 };
}

function calcSteam(T, P) {
  if (T < 273.15 || T > 1073.15) return { err:"T range: 273.15-1073.15 K" };
  if (P <= 0 || P > 100e6) return { err:"P range: 0-100 MPa" };
  // IAPWS-IF97 B23 boundary between Region 2 and Region 3:
  // P_B23(T) = n1 + n2*T + n3*T^2 (MPa), valid 623.15K-863.15K
  // Only warn when T>623.15K AND P is above B23 line (i.e. in Region 3)
  if (T > 623.15) {
    var n1=0.34805185628969e3, n2=-0.11671859879975e1, n3=0.10192970039326e-2;
    var P_B23 = (n1 + n2*T + n3*T*T) * 1e6; // Pa
    if (P >= P_B23) {
      return {err: "T="+((T-273.15).toFixed(1))+"\u00b0C, P="+(P/1e6).toFixed(2)+" MPa(a): IAPWS Region 3 (near critical). Not yet implemented."};
    }
  }
  var Ps = steamPsat(T);
  if (T <= 623.15 && Ps) {
    if (Math.abs(P-Ps)/Ps < 0.002) { var lq=steamR1(T,P), vp=steamR2(T,P); return { ph:"two-phase", phCN:"\u4e24\u76f8(\u9971\u548c)", Psat:Ps, liq:lq, vap:vp, dhv:vp.h-lq.h, rg:4 }; }
    if (P > Ps) { var r1=steamR1(T,P); r1.Psat=Ps; return r1; }
    var r2=steamR2(T,P); r2.Psat=Ps; return r2;
  }
  var r = steamR2(T, P);
  if (T > 647.096) { r.ph = "supercritical"; r.phCN = "\u8d85\u4e34\u754c"; }
  return r;
}

// ============================================================================
// PHASE 2: VLE - Vapor-Liquid Equilibrium for Multi-Component Mixtures
// Reference: Smith, Van Ness & Abbott, "Introduction to Chemical Engineering
//            Thermodynamics", 8th Ed, Ch. 10-14
// ============================================================================

// Mixture parameters for PR EOS (van der Waals one-fluid mixing rules)
// am = ΣΣ xi*xj*aij, where aij = sqrt(ai*aj)*(1-kij)
// bm = Σ xi*bi
function mixPR(T, x, clist) {
  var nc = clist.length;
  var ai = [], bi = [];
  for (var i = 0; i < nc; i++) {
    var ab = prAB(T, clist[i].Tc, clist[i].Pc, clist[i].w);
    ai.push(ab.a); bi.push(ab.b);
  }
  var am = 0, bm = 0;
  var aij = [];
  for (var ii = 0; ii < nc; ii++) {
    aij.push([]);
    for (var jj = 0; jj < nc; jj++) {
      var kij = getKij(clist[ii].id, clist[jj].id);
      aij[ii].push(Math.sqrt(ai[ii] * ai[jj]) * (1 - kij));
      am += x[ii] * x[jj] * aij[ii][jj];
    }
    bm += x[ii] * bi[ii];
  }
  return { am: am, bm: bm, ai: ai, bi: bi, aij: aij };
}

// Fugacity coefficient of component i in mixture (PR EOS)
// ln(φi) = bi/bm*(Z-1) - ln(Z-B) - A/(2√2*B) * [2*Σ(xj*aij)/am - bi/bm] * ln[...]
function lnPhiMix_PR(i, T, P, x, Z, clist) {
  var mx = mixPR(T, x, clist);
  var A = mx.am * P / (R_GAS * R_GAS * T * T);
  var B = mx.bm * P / (R_GAS * T);
  var s2 = 1.4142135623730951;
  var sumAij = 0;
  for (var j = 0; j < clist.length; j++) { sumAij += x[j] * mx.aij[i][j]; }
  var n = Z + (1 + s2) * B, d = Z + (1 - s2) * B;
  if (n <= 0 || d <= 0 || Z - B <= 0) return 30; // large positive = very low fugacity
  var lnTerm = Math.log(n / d);
  return mx.bi[i] / mx.bm * (Z - 1) - Math.log(Z - B) - A / (2 * s2 * B) * (2 * sumAij / mx.am - mx.bi[i] / mx.bm) * lnTerm;
}

// Solve cubic for mixture Z
function calcZmix(T, P, x, clist, phase) {
  var mx = mixPR(T, x, clist);
  var A = mx.am * P / (R_GAS * R_GAS * T * T);
  var B = mx.bm * P / (R_GAS * T);
  var roots = prRoots(A, B);
  if (roots.length === 0) return 0.9; // fallback
  return phase === "liquid" ? roots[0] : roots[roots.length - 1];
}

// Wilson K-value initial estimate
// Ki = (Pci/P) * exp[5.373*(1+ωi)*(1-Tci/T)]
function wilsonK(T, P, clist) {
  var K = [];
  for (var i = 0; i < clist.length; i++) {
    K.push(clist[i].Pc / P * Math.exp(5.373 * (1 + clist[i].w) * (1 - clist[i].Tc / T)));
  }
  return K;
}

// Update K from fugacity coefficients: Ki = φiL / φiV
function updateK(T, P, x, y, clist) {
  var nc = clist.length;
  var Zl = calcZmix(T, P, x, clist, "liquid");
  var Zv = calcZmix(T, P, y, clist, "vapor");
  var K = [];
  for (var i = 0; i < nc; i++) {
    var lpL = lnPhiMix_PR(i, T, P, x, Zl, clist);
    var lpV = lnPhiMix_PR(i, T, P, y, Zv, clist);
    K.push(Math.exp(lpL - lpV));
  }
  return K;
}

// Rachford-Rice equation: Σ zi*(Ki-1)/(1+β*(Ki-1)) = 0
function rachfordRice(beta, z, K) {
  var sum = 0;
  for (var i = 0; i < z.length; i++) {
    sum += z[i] * (K[i] - 1) / (1 + beta * (K[i] - 1));
  }
  return sum;
}

// TP Flash: given T, P, z → find β (vapor fraction), x, y
function tpFlash(T, P, z, clist, maxIter) {
  maxIter = maxIter || 100;
  var nc = z.length;
  var K = wilsonK(T, P, clist);

  // Check if single phase
  var sumKz = 0, sumZK = 0;
  for (var i = 0; i < nc; i++) { sumKz += K[i] * z[i]; sumZK += z[i] / K[i]; }
  if (sumKz <= 1) return { beta: 0, x: z.slice(), y: z.map(function(zi, i) { return K[i] * zi; }), K: K, ph: "liquid", iter: 0 };
  if (sumZK <= 1) return { beta: 1, x: z.map(function(zi, i) { return zi / K[i]; }), y: z.slice(), K: K, ph: "vapor", iter: 0 };

  // Solve Rachford-Rice for beta using Brent method
  var beta = 0.5;
  for (var outer = 0; outer < maxIter; outer++) {
    // Find beta bounds
    var bLo = 0, bHi = 1;
    // Tighten bounds to avoid division by zero
    for (var ib = 0; ib < nc; ib++) {
      var val = 1 / (1 - K[ib]);
      if (K[ib] > 1 && val > bLo) bLo = Math.max(bLo, val + 1e-10);
      if (K[ib] < 1 && val < bHi) bHi = Math.min(bHi, val - 1e-10);
    }
    bLo = Math.max(bLo, 1e-10); bHi = Math.min(bHi, 1 - 1e-10);
    if (bLo >= bHi) bLo = 1e-10;

    // Bisection for beta
    for (var bis = 0; bis < 80; bis++) {
      beta = (bLo + bHi) / 2;
      var rr = rachfordRice(beta, z, K);
      if (Math.abs(rr) < 1e-12) break;
      if (rr > 0) bLo = beta; else bHi = beta;
    }

    // Compute x, y from beta and K
    var x = [], y = [];
    for (var ix = 0; ix < nc; ix++) {
      x.push(z[ix] / (1 + beta * (K[ix] - 1)));
      y.push(K[ix] * z[ix] / (1 + beta * (K[ix] - 1)));
    }

    // Normalize
    var sx = 0, sy = 0;
    for (var in2 = 0; in2 < nc; in2++) { sx += x[in2]; sy += y[in2]; }
    for (var in3 = 0; in3 < nc; in3++) { x[in3] /= sx; y[in3] /= sy; }

    // Update K from EOS
    var Knew = updateK(T, P, x, y, clist);
    var maxDiff = 0;
    for (var ik = 0; ik < nc; ik++) {
      if (K[ik] > 0) maxDiff = Math.max(maxDiff, Math.abs(Knew[ik] / K[ik] - 1));
      K[ik] = Knew[ik];
    }
    if (maxDiff < 1e-8) {
      return { beta: beta, x: x, y: y, K: K, ph: "two-phase", iter: outer + 1 };
    }

    // Recheck phase
    sumKz = 0; sumZK = 0;
    for (var ic = 0; ic < nc; ic++) { sumKz += K[ic] * z[ic]; sumZK += z[ic] / K[ic]; }
    if (sumKz <= 1) return { beta: 0, x: z.slice(), y: z.map(function(zi2, i2) { return K[i2] * zi2; }), K: K, ph: "liquid", iter: outer + 1 };
    if (sumZK <= 1) return { beta: 1, x: z.map(function(zi3, i3) { return zi3 / K[i3]; }), y: z.slice(), K: K, ph: "vapor", iter: outer + 1 };
  }
  return { beta: beta, x: x || z, y: y || z, K: K, ph: "two-phase", iter: maxIter, warn: "max iter" };
}

// Bubble Point Pressure: given T, x → find P, y (bisection on Wilson + EOS refinement)
function bubbleP(T, x, clist) {
  var nc = x.length;
  // Bisection: find P such that Σ xi*Ki(T,P) = 1
  var Plo = 1000, Phi = 0;
  for (var i = 0; i < nc; i++) {
    var Psi = clist[i].Pc * Math.exp(5.373 * (1 + clist[i].w) * (1 - clist[i].Tc / T));
    if (Psi > Phi) Phi = Psi;
    if (Psi < Plo && Psi > 0) Plo = Psi;
  }
  Plo = Math.max(Plo * 0.01, 100); Phi = Math.min(Phi * 5, 1e9);
  var P = (Plo + Phi) / 2;
  for (var bis = 0; bis < 100; bis++) {
    P = (Plo + Phi) / 2;
    var K = wilsonK(T, P, clist);
    var s = 0; for (var j = 0; j < nc; j++) s += x[j] * K[j];
    if (Math.abs(s - 1) < 1e-10) break;
    if (s > 1) Plo = P; else Phi = P;
  }
  var Kw = wilsonK(T, P, clist);
  var y = []; var sy = 0;
  for (var jw = 0; jw < nc; jw++) { y.push(x[jw] * Kw[jw]); sy += y[jw]; }
  for (var jn = 0; jn < nc; jn++) y[jn] /= sy;
  // EOS refinement
  for (var iter = 0; iter < 50; iter++) {
    var Zl2 = calcZmix(T, P, x, clist, "liquid");
    var Zv2 = calcZmix(T, P, y, clist, "vapor");
    var Ke = [], s2 = 0;
    for (var ie2 = 0; ie2 < nc; ie2++) { Ke.push(Math.exp(lnPhiMix_PR(ie2,T,P,x,Zl2,clist)-lnPhiMix_PR(ie2,T,P,y,Zv2,clist))); }
    y = []; s2 = 0;
    for (var je = 0; je < nc; je++) { y.push(x[je]*Ke[je]); s2+=y[je]; }
    for (var je2=0;je2<nc;je2++) y[je2]/=s2;
    var Pnew = P * s2;
    if (Math.abs(Pnew-P)/P < 1e-8) return { P:Pnew, y:y, K:Ke, iter:iter+1 };
    P = P + 0.5*(Pnew-P); if (P<=0) P=100;
  }
  return { P:P, y:y, K:Kw, iter:50 };
}

// Bubble Point Temperature: given P, x → find T, y (bisection + EOS refinement)
function bubbleT(P, x, clist) {
  var nc = x.length;
  var Tlo = 10, Thi = 0;
  for (var i = 0; i < nc; i++) {
    var Tbi = clist[i].Tc / (1 - Math.log(P / clist[i].Pc) / (5.373 * (1 + clist[i].w)));
    if (isNaN(Tbi) || Tbi <= 0) Tbi = clist[i].Tc * 0.7;
    if (Tbi > Thi) Thi = Tbi; if (Tbi < Tlo) Tlo = Tbi;
  }
  Tlo = Math.max(Tlo * 0.8, 10); Thi = Math.min(Thi * 1.2, 2000);
  var T = (Tlo+Thi)/2;
  // Bisection with Wilson
  for (var bis = 0; bis < 100; bis++) {
    T = (Tlo+Thi)/2;
    var K=wilsonK(T,P,clist); var s=0;
    for (var j=0;j<nc;j++) s+=x[j]*K[j];
    if (Math.abs(s-1)<1e-10) break;
    if (s>1) Thi=T; else Tlo=T;
  }
  var Kw=wilsonK(T,P,clist); var y=[]; var sy=0;
  for (var jw=0;jw<nc;jw++){y.push(x[jw]*Kw[jw]);sy+=y[jw]}
  for (var jn=0;jn<nc;jn++) y[jn]/=sy;
  var T_wilson=T, y_wilson=y.slice();
  // EOS refinement with divergence guard
  for (var iter=0;iter<50;iter++){
    var Zl=calcZmix(T,P,x,clist,"liquid");var Zv=calcZmix(T,P,y,clist,"vapor");
    var Ke=[],s2=0,ok=true;
    for(var ie=0;ie<nc;ie++){
      var lp=lnPhiMix_PR(ie,T,P,x,Zl,clist)-lnPhiMix_PR(ie,T,P,y,Zv,clist);
      if(!isFinite(lp)){ok=false;break;}
      Ke.push(Math.exp(Math.max(-30,Math.min(30,lp))));
    }
    if(!ok) return{T:T_wilson,y:y_wilson,K:Kw,iter:iter+1};
    y=[];s2=0;for(var je=0;je<nc;je++){y.push(x[je]*Ke[je]);s2+=y[je]}
    if(s2<=0||!isFinite(s2)) return{T:T_wilson,y:y_wilson,K:Kw,iter:iter+1};
    for(var je2=0;je2<nc;je2++) y[je2]/=s2;
    if(Math.abs(s2-1)<1e-8) return{T:T,y:y,K:Ke,iter:iter+1};
    var dT=0.01;
    var Zl2=calcZmix(T+dT,P,x,clist,"liquid");var Zv2=calcZmix(T+dT,P,y,clist,"vapor");
    var s3=0;
    for(var je3=0;je3<nc;je3++){
      var lp3=lnPhiMix_PR(je3,T+dT,P,x,Zl2,clist)-lnPhiMix_PR(je3,T+dT,P,y,Zv2,clist);
      s3+=x[je3]*Math.exp(Math.max(-30,Math.min(30,lp3)));
    }
    var deriv=(s3-s2)/dT;
    if(Math.abs(deriv)>1e-12){
      var step=-(s2-1)/deriv;
      if(Math.abs(step)>5)step=5*(step>0?1:-1);
      var Tnew=T+0.7*step;
      if(Tnew<10||Tnew>2000||!isFinite(Tnew)) return{T:T_wilson,y:y_wilson,K:Kw,iter:iter+1};
      T=Tnew;
    }
  }
  return{T:T,y:y,K:Kw,iter:50};
}

// Dew Point Pressure: given T, y → find P, x
function dewP(T, y, clist) {
  var nc = y.length;
  var Plo = 1000, Phi = 0;
  for (var i = 0; i < nc; i++) {
    var Psi = clist[i].Pc * Math.exp(5.373*(1+clist[i].w)*(1-clist[i].Tc/T));
    if (Psi > Phi) Phi = Psi; if (Psi < Plo && Psi > 0) Plo = Psi;
  }
  Plo = Math.max(Plo*0.01, 100); Phi = Math.min(Phi*5, 1e9);
  var P=(Plo+Phi)/2;
  for (var bis=0;bis<100;bis++){
    P=(Plo+Phi)/2;var K=wilsonK(T,P,clist);var s=0;for(var j=0;j<nc;j++)s+=y[j]/K[j];
    if(Math.abs(s-1)<1e-10)break;
    if(s<1)Plo=P;else Phi=P;
  }
  var Kw=wilsonK(T,P,clist);var x=[];var sx=0;
  for(var jw=0;jw<nc;jw++){x.push(y[jw]/Kw[jw]);sx+=x[jw]}
  for(var jn=0;jn<nc;jn++)x[jn]/=sx;
  for(var iter=0;iter<50;iter++){
    var Zl=calcZmix(T,P,x,clist,"liquid");var Zv=calcZmix(T,P,y,clist,"vapor");
    var Ke=[];for(var ie=0;ie<nc;ie++)Ke.push(Math.exp(lnPhiMix_PR(ie,T,P,x,Zl,clist)-lnPhiMix_PR(ie,T,P,y,Zv,clist)));
    var s2=0;x=[];for(var je=0;je<nc;je++){x.push(y[je]/Ke[je]);s2+=x[je]}
    for(var je2=0;je2<nc;je2++)x[je2]/=s2;
    var Pnew=P/s2;
    if(Math.abs(Pnew-P)/P<1e-8)return{P:Pnew,x:x,K:Ke,iter:iter+1};
    P=P+0.5*(Pnew-P);if(P<=0)P=100;
  }
  return{P:P,x:x,K:Kw,iter:50};
}

// Dew Point Temperature: given P, y → find T, x (bisection + EOS)
function dewT(P, y, clist) {
  var nc = y.length;
  var Tlo=10,Thi=0;
  for(var i=0;i<nc;i++){var Tbi=clist[i].Tc/(1-Math.log(P/clist[i].Pc)/(5.373*(1+clist[i].w)));if(isNaN(Tbi)||Tbi<=0)Tbi=clist[i].Tc*0.7;if(Tbi>Thi)Thi=Tbi;if(Tbi<Tlo)Tlo=Tbi}
  Tlo=Math.max(Tlo*0.8,10);Thi=Math.min(Thi*1.2,2000);var T=(Tlo+Thi)/2;
  for(var bis=0;bis<100;bis++){
    T=(Tlo+Thi)/2;var K=wilsonK(T,P,clist);var s=0;for(var j=0;j<nc;j++)s+=y[j]/K[j];
    if(Math.abs(s-1)<1e-10)break;
    if(s<1)Thi=T;else Tlo=T;
  }
  var Kw=wilsonK(T,P,clist);var x=[];var sx=0;
  for(var jw=0;jw<nc;jw++){x.push(y[jw]/Kw[jw]);sx+=x[jw]}
  for(var jn=0;jn<nc;jn++)x[jn]/=sx;
  var T_wilson=T, x_wilson=x.slice();
  // EOS refinement with divergence guard
  for(var iter=0;iter<50;iter++){
    var Zl=calcZmix(T,P,x,clist,"liquid");var Zv=calcZmix(T,P,y,clist,"vapor");
    var Ke=[],ok=true;
    for(var ie=0;ie<nc;ie++){
      var lp=lnPhiMix_PR(ie,T,P,x,Zl,clist)-lnPhiMix_PR(ie,T,P,y,Zv,clist);
      if(!isFinite(lp)){ok=false;break;}
      Ke.push(Math.exp(Math.max(-30,Math.min(30,lp))));
    }
    if(!ok) return{T:T_wilson,x:x_wilson,K:Kw,iter:iter+1};
    var s2=0;x=[];for(var je=0;je<nc;je++){x.push(y[je]/Ke[je]);s2+=x[je]}
    if(s2<=0||!isFinite(s2)) return{T:T_wilson,x:x_wilson,K:Kw,iter:iter+1};
    for(var je2=0;je2<nc;je2++)x[je2]/=s2;
    if(Math.abs(s2-1)<1e-8)return{T:T,x:x,K:Ke,iter:iter+1};
    var dT=0.01;
    var Zl2=calcZmix(T+dT,P,x,clist,"liquid");var Zv2=calcZmix(T+dT,P,y,clist,"vapor");
    var s3=0;
    for(var je3=0;je3<nc;je3++){
      var lp3=lnPhiMix_PR(je3,T+dT,P,x,Zl2,clist)-lnPhiMix_PR(je3,T+dT,P,y,Zv2,clist);
      s3+=y[je3]/Math.exp(Math.max(-30,Math.min(30,lp3)));
    }
    var deriv=(s3-s2)/dT;
    if(Math.abs(deriv)>1e-12){
      var step=-(s2-1)/deriv;
      if(Math.abs(step)>5)step=5*(step>0?1:-1);
      var Tnew=T+0.7*step;
      if(Tnew<10||Tnew>2000||!isFinite(Tnew)) return{T:T_wilson,x:x_wilson,K:Kw,iter:iter+1};
      T=Tnew;
    }
  }
  return{T:T,x:x,K:wilsonK(T,P,clist),iter:50};
}

// Mixture H/S departure via PR EOS
// phase: "liquid" -> smallest Z root; "vapor" -> largest Z root; default "vapor"
function mixHSdep(T, P, x, clist, phase) {
  var nc = clist.length;
  var mx = mixPR(T, x, clist);
  var A = mx.am * P / (R_GAS * R_GAS * T * T);
  var B = mx.bm * P / (R_GAS * T);
  var roots = prRoots(A, B);
  if (!roots || roots.length === 0) return {Hd:0, Sd:0};
  var Z = (phase === "liquid" && roots.length > 1) ? roots[0] : roots[roots.length - 1];
  var s2 = 1.4142135623730951;
  // Precompute per-component sqrt(ai) and dai/dT for efficiency
  var sqAi = [], daidT = [];
  for (var i = 0; i < nc; i++) {
    var abi = prAB(T, clist[i].Tc, clist[i].Pc, clist[i].w);
    var mi = 0.37464 + 1.54226*clist[i].w - 0.26992*clist[i].w*clist[i].w;
    var Tri = T / clist[i].Tc;
    sqAi.push(Math.sqrt(abi.a));
    // da_i/dT = -0.45724*R²Tc²/Pc * m_i * sqrt(alpha_i) / (sqrt(Tr_i)*Tc_i)
    // sqrt(alpha_i) = sqrt(abi.a) / sqrt(0.45724*R²Tc²/Pc) = sqAi[i]/sqrt(ac_i)
    var ac_i = 0.45724 * R_GAS * R_GAS * clist[i].Tc * clist[i].Tc / clist[i].Pc;
    daidT.push(-ac_i * mi * sqAi[i] / (Math.sqrt(ac_i) * Math.sqrt(Tri) * clist[i].Tc));
  }
  // dam/dT = ΣΣ xi*xj*(1-kij)*d(sqrt(ai*aj))/dT
  //        = ΣΣ xi*xj*(1-kij)*(daidT*sqAj + dajdT*sqAi)/(2*sqAi*sqAj) * sqAi*sqAj ... simplifies to:
  // dam/dT = ΣΣ xi*xj*(1-kij)*(daidT[i]*sqAi[j] + daidT[j]*sqAi[i]) / 2
  var dadT = 0;
  for (var ii = 0; ii < nc; ii++) {
    for (var jj = 0; jj < nc; jj++) {
      var kij = getKij(clist[ii].id, clist[jj].id);
      dadT += x[ii] * x[jj] * (1 - kij) * (daidT[ii] * sqAi[jj] + daidT[jj] * sqAi[ii]) / 2;
    }
  }
  var n = Z + (1 + s2) * B, d = Z + (1 - s2) * B;
  if (n <= 0 || d <= 0 || Z - B <= 0) return {Hd:0, Sd:0};
  var lnT = Math.log(n / d);
  var Hd = R_GAS * T * (Z - 1) + (T * dadT - mx.am) / (2 * s2 * mx.bm) * lnT;
  var Sd = R_GAS * Math.log(Math.max(Z - B, 1e-30)) + dadT / (2 * s2 * mx.bm) * lnT;
  return {Hd: Hd, Sd: Sd, Z: Z};
}

// Total mixture molar enthalpy (J/mol)
function mixEnthalpy(T, P, x, clist, phase) {
  var Hig = 0;
  for (var i = 0; i < clist.length; i++) Hig += x[i] * hIg(T, clist[i]);
  var dep = mixHSdep(T, P, x, clist, phase || "vapor");
  return Hig + dep.Hd;
}

// Total mixture molar entropy (J/(mol*K))
function mixEntropy(T, P, x, clist, phase) {
  var Sig = 0;
  for (var i = 0; i < clist.length; i++) Sig += x[i] * sIg(T, P, clist[i]);
  var dep = mixHSdep(T, P, x, clist, phase || "vapor");
  return Sig + dep.Sd;
}

// VLE PH Flash: given P(Pa), H_target(J/mol), z, clist → find T via outer Newton
// Uses analytical dH/dT ≈ Σzi*Cp,ig,i (ideal gas heat capacity) to avoid double tpFlash call
function vlePhFlash(P, H_target, z, clist) {
  // Better initial T: estimate from ideal gas H = Σzi*hIg(T,ci) ≈ H_target
  // Use Cp at 300K to get rough T0
  var Cp0 = 0;
  for (var ic = 0; ic < clist.length; ic++) Cp0 += z[ic] * R_GAS * (clist[ic].cp[0] + clist[ic].cp[1]*300 + clist[ic].cp[2]*300*300);
  var H300 = 0;
  for (var ic2 = 0; ic2 < clist.length; ic2++) H300 += z[ic2] * hIg(300, clist[ic2]);
  var T = Math.max(120, Math.min(900, 300 + (H_target - H300) / Math.max(Cp0, 10)));
  for (var oit = 0; oit < 60; oit++) {
    var fl = tpFlash(T, P, z, clist);
    var beta = fl.beta !== undefined ? fl.beta : 0;
    var xl = fl.x || z, xv = fl.y || z;
    // Correct phase for Z root selection
    var HL = mixEnthalpy(T, P, xl, clist, "liquid");
    var HV = mixEnthalpy(T, P, xv, clist, "vapor");
    var Hmix = (1 - beta) * HL + beta * HV;
    var diff = Hmix - H_target;
    if (Math.abs(diff) < 0.5) return Object.assign({T: T, Hmix: Hmix}, fl);
    // Analytical dH/dT ≈ mixture Cp (ideal gas + small departure correction)
    var dHdT = 0;
    for (var i = 0; i < clist.length; i++) dHdT += z[i] * R_GAS * (clist[i].cp[0] + clist[i].cp[1]*T + clist[i].cp[2]*T*T + clist[i].cp[3]*T*T*T);
    dHdT = Math.max(dHdT, 5); // floor to prevent zero
    var step = -diff / dHdT;
    if (Math.abs(step) > 40) step = 40 * (step > 0 ? 1 : -1);
    T += step;
    if (T < 100) T = 100;
    if (T > 1000) T = 1000;
  }
  var flFinal = tpFlash(T, P, z, clist);
  return Object.assign({T: T, warn: "max iter"}, flFinal);
}

// VLE PS Flash: given P(Pa), S_target(J/(mol*K)), z, clist → find T
function vlePsFlash(P, S_target, z, clist) {
  // Initial estimate from ideal gas entropy: S_ig(T,P) ≈ S_target + entropy correction
  var T = 300;
  for (var oit = 0; oit < 60; oit++) {
    var fl = tpFlash(T, P, z, clist);
    var beta = fl.beta !== undefined ? fl.beta : 0;
    var xl = fl.x || z, xv = fl.y || z;
    var SL = mixEntropy(T, P, xl, clist, "liquid");
    var SV = mixEntropy(T, P, xv, clist, "vapor");
    var Smix = (1 - beta) * SL + beta * SV;
    var diff = Smix - S_target;
    if (Math.abs(diff) < 0.001) return Object.assign({T: T, Smix: Smix}, fl);
    // dS/dT ≈ Cp/T (ideal gas)
    var Cp = 0;
    for (var i = 0; i < clist.length; i++) Cp += z[i] * R_GAS * (clist[i].cp[0] + clist[i].cp[1]*T + clist[i].cp[2]*T*T + clist[i].cp[3]*T*T*T);
    var dSdT = Math.max(Cp / T, 0.005);
    var step = -diff / dSdT;
    if (Math.abs(step) > 40) step = 40 * (step > 0 ? 1 : -1);
    T += step;
    if (T < 100) T = 100;
    if (T > 1000) T = 1000;
  }
  var flFinal = tpFlash(T, P, z, clist);
  return Object.assign({T: T, warn: "max iter"}, flFinal);
}


function genTxy(P, c1, c2, nPts) {
  nPts = nPts || 30;
  var clist = [c1, c2];
  // Physical range: estimate Tsat for each pure component at given P (Lee-Kesler)
  var TsatEst = function(comp) {
    var Tr_lo=0.35, Tr_hi=0.999, P_lo=comp.Pc*Math.exp(5.92714-6.09648/Tr_lo-1.28862*Math.log(Tr_lo)+0.169347*Math.pow(Tr_lo,6)+comp.w*(15.2518-15.6875/Tr_lo-13.4721*Math.log(Tr_lo)+0.43577*Math.pow(Tr_lo,6)));
    if (P < P_lo) return comp.Tc * 0.35; // below triple-point pressure estimate
    for (var bi=0; bi<80; bi++) {
      var Tr=(Tr_lo+Tr_hi)/2;
      var f0=5.92714-6.09648/Tr-1.28862*Math.log(Tr)+0.169347*Math.pow(Tr,6);
      var f1=15.2518-15.6875/Tr-13.4721*Math.log(Tr)+0.43577*Math.pow(Tr,6);
      var Ps=comp.Pc*Math.exp(f0+comp.w*f1);
      if(Ps<P) Tr_lo=Tr; else Tr_hi=Tr;
      if(Math.abs(Ps-P)/P<1e-6) break;
    }
    return ((Tr_lo+Tr_hi)/2)*comp.Tc;
  };
  var Ts1=TsatEst(c1), Ts2=TsatEst(c2);
  var Tmin=Math.min(Ts1,Ts2)*0.95 - 5; // K, with margin
  var Tmax=Math.max(Ts1,Ts2)*1.05 + 5;
  // Warn if system is likely problematic
  var supercrit = (Ts1 < 0.4*c1.Tc) || (Ts2 < 0.4*c2.Tc);
  var data = [];
  for (var i = 0; i <= nPts; i++) {
    var x1 = i / nPts;
    var bub = bubbleT(P, [x1, 1-x1], clist);
    var Tb = bub.T && !isNaN(bub.T) ? bub.T : null;
    // Filter physically unreasonable values
    if (Tb && (Tb < Tmin || Tb > Tmax)) Tb = null;
    var Tc_val = Tb ? Math.round((Tb-273.15)*100)/100 : null;
    var y1 = (bub.y && Tb) ? Math.round(bub.y[0]*1000)/1000 : x1;
    data.push({ x1: Math.round(x1*1000)/1000, Tbub: Tc_val, y1: y1 });
  }
  var dewData = [];
  for (var j = 0; j <= nPts; j++) {
    var yy1 = j / nPts;
    var dew = dewT(P, [yy1, 1-yy1], clist);
    var Td = dew.T && !isNaN(dew.T) ? dew.T : null;
    if (Td && (Td < Tmin || Td > Tmax)) Td = null;
    var Td_val = Td ? Math.round((Td-273.15)*100)/100 : null;
    dewData.push({ y1: Math.round(yy1*1000)/1000, Tdew: Td_val });
  }
  return { bub: data, dew: dewData, supercrit: supercrit,
    Ts1: Math.round((Ts1-273.15)*10)/10, Ts2: Math.round((Ts2-273.15)*10)/10 };
}

// Generate P-x-y diagram
function genPxy(T, c1, c2, nPts) {
  nPts = nPts || 30;
  var clist = [c1, c2];
  // Physical range: Psat of each pure component at T (Lee-Kesler)
  var PsatEst = function(comp) {
    var Tr = T / comp.Tc; if (Tr >= 1) return null;
    var f0=5.92714-6.09648/Tr-1.28862*Math.log(Tr)+0.169347*Math.pow(Tr,6);
    var f1=15.2518-15.6875/Tr-13.4721*Math.log(Tr)+0.43577*Math.pow(Tr,6);
    return comp.Pc * Math.exp(f0 + comp.w * f1);
  };
  var Ps1=PsatEst(c1), Ps2=PsatEst(c2);
  // Supercritical if one component has T > Tc
  var supercrit = (T >= c1.Tc) || (T >= c2.Tc);
  var Pmin=0, Pmax=Math.max(Ps1||0,Ps2||0)*1.3 + 1e5;
  var bubData = [], dewData = [];
  for (var i = 0; i <= nPts; i++) {
    var x1 = i / nPts;
    var bub = bubbleP(T, [x1, 1-x1], clist);
    var Pb = bub.P && !isNaN(bub.P) && bub.P > 0 && bub.P < Pmax ? Math.round(bub.P/1e6*10000)/10000 : null;
    bubData.push({ x1: Math.round(x1*1000)/1000, Pbub: Pb, y1: bub.y ? Math.round(bub.y[0]*1000)/1000 : x1 });
  }
  for (var j = 0; j <= nPts; j++) {
    var yy1 = j / nPts;
    var dew = dewP(T, [yy1, 1-yy1], clist);
    var Pd = dew.P && !isNaN(dew.P) && dew.P > 0 && dew.P < Pmax ? Math.round(dew.P/1e6*10000)/10000 : null;
    dewData.push({ y1: Math.round(yy1*1000)/1000, Pdew: Pd });
  }
  return { bub: bubData, dew: dewData, supercrit: supercrit,
    Ps1: Ps1 ? Math.round(Ps1/1e6*1000)/1000 : null,
    Ps2: Ps2 ? Math.round(Ps2/1e6*1000)/1000 : null };
}


// ============================================================================
// EQUIPMENT CALCULATIONS: Compressor, Pump, Pipe
// ============================================================================

// Compressor: Isentropic compression with efficiency correction
// Input: P1(Pa), T1(K), P2(Pa), eta_s (isentropic eff), gas Cp/Cv(gamma), MW
// Output: T2_ideal, T2_actual, Work, Power
// Method: T2s = T1*(P2/P1)^((gamma-1)/gamma), T2 = T1 + (T2s-T1)/eta
function calcCompressor(P1, T1, P2, etaS, gamma, MW, massFlow) {
  if (P2 <= P1) return {err:"P2 must > P1"};
  if (etaS <= 0 || etaS > 1) return {err:"Efficiency 0-1"};
  var ratio = P2/P1;
  var exp1 = (gamma-1)/gamma;
  var T2s = T1 * Math.pow(ratio, exp1); // isentropic outlet T
  var T2 = T1 + (T2s - T1) / etaS; // actual outlet T
  var Cp = gamma * R_GAS / ((gamma-1) * MW/1000); // J/(kg*K)
  var ws = Cp * (T2s - T1); // isentropic specific work J/kg
  var wa = Cp * (T2 - T1); // actual specific work J/kg
  var power = massFlow > 0 ? wa * massFlow / 1000 : 0; // kW
  return {T2s:T2s, T2:T2, ws:ws/1000, wa:wa/1000, power:power, ratio:ratio, polyExp:Math.log(ratio)/Math.log(T2/T1)*exp1};
}

// Expander/Turbine: Isentropic expansion
// P1 > P2 (inlet pressure > outlet pressure), eta = isentropic efficiency
function calcExpander(P1, T1, P2, etaS, gamma, MW, massFlow) {
  if (P1 <= P2) return {err:"P1 must > P2 (expansion)"};
  if (etaS <= 0 || etaS > 1) return {err:"Efficiency 0-1"};
  var ratio = P1/P2; // expansion ratio > 1
  var exp1 = (gamma-1)/gamma;
  var T2s = T1 / Math.pow(ratio, exp1); // isentropic outlet T (lower)
  var T2 = T1 - (T1 - T2s) * etaS; // actual outlet T (less drop due to irreversibility)
  var Cp = gamma * R_GAS / ((gamma-1) * MW/1000); // J/(kg*K)
  var ws = Cp * (T1 - T2s); // isentropic specific work output J/kg (positive)
  var wa = Cp * (T1 - T2); // actual specific work output J/kg
  var power = massFlow > 0 ? wa * massFlow / 1000 : 0; // kW output
  return {T2s:T2s, T2:T2, ws:ws/1000, wa:wa/1000, power:power, ratio:ratio};
}

// Pump: Hydraulic power + efficiency
// Input: flow(m3/h), dP(Pa), eta(pump eff), rho(kg/m3)
function calcPump(flow_m3h, P1, P2, eta, rho) {
  if (P2 <= P1) return {err:"P2 must > P1"};
  var dP = P2 - P1;
  var Q = flow_m3h / 3600; // m3/s
  var head = dP / (rho * 9.80665); // m
  var Ph = rho * 9.80665 * head * Q; // W hydraulic
  var Ps = Ph / eta; // W shaft
  return {dP:dP, head:head, Ph:Ph/1000, Ps:Ps/1000, Q:Q};
}

// Pipe pressure drop: Darcy-Weisbach
// Input: L(m), D(m), roughness(m), velocity(m/s) or flow, rho, mu(Pa.s)
function calcPipe(L, D, eps, vel, rho, mu) {
  var Re = rho * vel * D / mu;
  var f;
  if (Re < 2300) {
    f = 64 / Re; // laminar
  } else {
    // Colebrook-White (explicit Swamee-Jain approximation)
    var A = eps / (3.7 * D);
    var B = 5.74 / Math.pow(Re, 0.9);
    f = 0.25 / Math.pow(Math.log10(A + B), 2);
  }
  var dP = f * L / D * rho * vel * vel / 2; // Pa
  var dH = dP / (rho * 9.80665); // m head loss
  return {Re:Re, f:f, dP:dP, dH:dH, regimeEn:Re<2300?"Laminar":"Turbulent",regime:Re<2300?"层流/Laminar":"湍流/Turbulent"};
}

// Unit conversion tables
// ============================================================================
// PH and PS FLASH for IAPWS-IF97 (backward equations)
// Given P+H or P+S, find T by Newton iteration on Region 2
// ============================================================================
function steamPHflash(P, H_target) {
  if (P <= 0 || P > 100e6) return {err:"P range: 0-100 MPa"};
  // B23 Region 3 check (same as calcSteam)
  if (P > 16.529e6) { var Tsat_b23 = iapws97_Tsat(P); if (!Tsat_b23 || Tsat_b23 > 623.15) return {err:"P="+( P/1e6).toFixed(2)+" MPa(a): IAPWS Region 3 not implemented. Use P<16.5 MPa(a)."}; }
  var Tsat = iapws97_Tsat(P);
  var hf = null, hg = null;
  if (Tsat && Tsat >= 273.15 && Tsat <= 623.15) {
    var satLiq = steamR1(Tsat, P);
    var satVap = steamR2(Tsat, P);
    hf = satLiq.h; hg = satVap.h;
    if (H_target >= hf && H_target <= hg) {
      var x = (H_target - hf)/(hg - hf);
      var sf=satLiq.s, sg=satVap.s, vf=satLiq.v, vg=satVap.v;
      return {T:Tsat, props:{ph:"two-phase",phCN:"\u4e24\u76f8(\u6e7f\u84b8\u6c7d)",h:H_target,s:sf+x*(sg-sf),v:vf+x*(vg-vf),rho:1/(vf+x*(vg-vf)),cp:NaN,cv:NaN,w:NaN,u:H_target-P*(vf+x*(vg-vf))/1000,rg:4,vaporFrac:x,Psat:P,Tsat:Tsat},iter:0};
    }
    if (H_target < hf) {
      var T = Tsat - 10;
      for (var iter=0; iter<100; iter++) {
        if (T<273.15) T=274; if (T>Tsat) T=Tsat-0.1;
        var r1 = steamR1(T, P);
        var diff = r1.h - H_target;
        if (Math.abs(diff) < 0.001) return {T:T,props:r1,iter:iter+1};
        var r1b = steamR1(T+0.1, P);
        var dhdT = (r1b.h-r1.h)/0.1;
        if (Math.abs(dhdT)<1e-10) dhdT=4.18;
        var step = -diff/dhdT;
        if (Math.abs(step)>20) step=20*(step>0?1:-1);
        T += step;
      }
      return {T:T,props:steamR1(T,P),iter:100,warn:"max iter"};
    }
  }
  var T2 = Tsat ? Tsat+20 : 500;
  for (var iter2=0; iter2<100; iter2++) {
    if (T2<273.15) T2=280; if (T2>1073) T2=1000;
    var r = steamR2(T2, P);
    if (!r||r.h===undefined) {T2+=10;continue;}
    var diff2 = r.h - H_target;
    if (Math.abs(diff2)<0.001) return {T:T2,props:r,iter:iter2+1};
    var r2b = steamR2(T2+0.1, P);
    var dhdT2 = (r2b.h-r.h)/0.1;
    if (Math.abs(dhdT2)<1e-10) dhdT2=2;
    var step2 = -diff2/dhdT2;
    if (Math.abs(step2)>50) step2=50*(step2>0?1:-1);
    T2 += step2;
  }
  return {T:T2,props:steamR2(T2,P),iter:100,warn:"max iter"};
}

function steamPSflash(P, S_target) {
  if (P <= 0 || P > 100e6) return {err:"P range: 0-100 MPa"};
  // B23 Region 3 check
  if (P > 16.529e6) { var Tsat_b23 = iapws97_Tsat(P); if (!Tsat_b23 || Tsat_b23 > 623.15) return {err:"P="+(P/1e6).toFixed(2)+" MPa(a): IAPWS Region 3 not implemented. Use P<16.5 MPa(a)."}; }
  var Tsat = iapws97_Tsat(P);
  var sf = null, sg = null;
  if (Tsat && Tsat >= 273.15 && Tsat <= 623.15) {
    var satLiq = steamR1(Tsat, P);
    var satVap = steamR2(Tsat, P);
    sf = satLiq.s; sg = satVap.s;
    if (S_target >= sf && S_target <= sg) {
      var x = (S_target-sf)/(sg-sf);
      var hf=satLiq.h, hg=satVap.h, vf=satLiq.v, vg=satVap.v;
      var h=hf+x*(hg-hf);
      return {T:Tsat,props:{ph:"two-phase",phCN:"\u4e24\u76f8(\u6e7f\u84b8\u6c7d)",h:h,s:S_target,v:vf+x*(vg-vf),rho:1/(vf+x*(vg-vf)),cp:NaN,cv:NaN,w:NaN,u:h-P*(vf+x*(vg-vf))/1000,rg:4,vaporFrac:x,Psat:P,Tsat:Tsat},iter:0};
    }
    if (S_target < sf) {
      var T = Tsat - 10;
      for (var iter=0; iter<100; iter++) {
        if (T<273.15) T=274; if (T>Tsat) T=Tsat-0.1;
        var r1 = steamR1(T, P);
        var diff = r1.s - S_target;
        if (Math.abs(diff)<0.0001) return {T:T,props:r1,iter:iter+1};
        var r1b = steamR1(T+0.1, P);
        var dsdT = (r1b.s-r1.s)/0.1;
        if (Math.abs(dsdT)<1e-10) dsdT=0.01;
        var step = -diff/dsdT;
        if (Math.abs(step)>20) step=20*(step>0?1:-1);
        T += step;
      }
      return {T:T,props:steamR1(T,P),iter:100,warn:"max iter"};
    }
  }
  var T2 = Tsat ? Tsat+20 : 500;
  for (var iter2=0; iter2<100; iter2++) {
    if (T2<273.15) T2=280; if (T2>1073) T2=1000;
    var r = steamR2(T2, P);
    if (!r||r.s===undefined) {T2+=10;continue;}
    var diff2 = r.s - S_target;
    if (Math.abs(diff2)<0.0001) return {T:T2,props:r,iter:iter2+1};
    var r2b = steamR2(T2+0.1, P);
    var dsdT2 = (r2b.s-r.s)/0.1;
    if (Math.abs(dsdT2)<1e-10) dsdT2=0.01;
    var step2 = -diff2/dsdT2;
    if (Math.abs(step2)>50) step2=50*(step2>0?1:-1);
    T2 += step2;
  }
  return {T:T2,props:steamR2(T2,P),iter:100,warn:"max iter"};
}

// ============================================================================
// REFRIGERATION CYCLE COP CALCULATOR
// Simple vapor compression cycle: evaporator → compressor → condenser → valve
// ============================================================================
function calcRefrigCOP(Tevap_K, Tcond_K, Tsh, Tsc, etaComp, comp) {
  // Tevap_K: evaporator saturation temp (K)
  // Tcond_K: condenser saturation temp (K)
  // Tsh: superheat at compressor inlet (K)
  // Tsc: subcooling at condenser outlet (K)
  // etaComp: compressor isentropic efficiency
  // comp: refrigerant component object
  // Carnot COP for reference
  var COPcarnot = Tevap_K / (Tcond_K - Tevap_K);
  // Estimate Pevap and Pcond using Lee-Kesler
  var Tr_e = Tevap_K / comp.Tc, Tr_c = Tcond_K / comp.Tc;
  var lkP = function(Tr) {
    var f0 = 5.92714 - 6.09648/Tr - 1.28862*Math.log(Tr) + 0.169347*Math.pow(Tr,6);
    var f1 = 15.2518 - 15.6875/Tr - 13.4721*Math.log(Tr) + 0.43577*Math.pow(Tr,6);
    return comp.Pc * Math.exp(f0 + comp.w * f1);
  };
  var Pevap = lkP(Tr_e), Pcond = lkP(Tr_c);
  if (!Pevap || !Pcond || Pevap <= 0 || Pcond <= 0) return {err:"Cannot estimate Psat"};
  var pressRatio = Pcond / Pevap;
  // Cp/Cv ratio from ideal gas Cp
  var CpIg = cpIg(Tevap_K + Tsh, comp); // J/(mol*K)
  var gamma = CpIg / (CpIg - R_GAS);
  // Compressor work (isentropic)
  var T1 = Tevap_K + Tsh; // compressor inlet
  var exp1 = (gamma - 1) / gamma;
  var T2s = T1 * Math.pow(pressRatio, exp1);
  var T2 = T1 + (T2s - T1) / etaComp;
  var Cp_mass = CpIg / (comp.MW / 1000); // J/(kg*K)
  var wComp = Cp_mass * (T2 - T1) / 1000; // kJ/kg
  // Evaporator capacity (approximate): qEvap ≈ latent heat
  // Use Pitzer correlation for Hvap: ΔHvap/RTc = 7.08*(1-Tr)^0.354 + 10.95*ω*(1-Tr)^0.456
  var dHvap_e = R_GAS * comp.Tc * (7.08*Math.pow(1-Tr_e,0.354) + 10.95*comp.w*Math.pow(1-Tr_e,0.456)); // J/mol
  var qEvap = dHvap_e / (comp.MW / 1000) / 1000; // kJ/kg
  var COP = qEvap / wComp;
  return {
    COP:COP, COPcarnot:COPcarnot, Pevap:Pevap, Pcond:Pcond, pressRatio:pressRatio,
    T1:T1, T2s:T2s, T2:T2, wComp:wComp, qEvap:qEvap, gamma:gamma
  };
}

// ============================================================================
// LMTD CALCULATOR (Log Mean Temperature Difference)
// ============================================================================
function calcLMTD(Th_in, Th_out, Tc_in, Tc_out, flowType) {
  // flowType: "counter" or "parallel"
  var dT1, dT2;
  if (flowType === "parallel") {
    dT1 = Th_in - Tc_in; // hot inlet - cold inlet
    dT2 = Th_out - Tc_out; // hot outlet - cold outlet
  } else {
    // counterflow (default)
    dT1 = Th_in - Tc_out; // hot inlet - cold outlet
    dT2 = Th_out - Tc_in; // hot outlet - cold inlet
  }
  if (dT1 <= 0 || dT2 <= 0) return {err:"Temperature cross detected"};
  if (Math.abs(dT1 - dT2) < 0.01) return {LMTD: dT1, dT1:dT1, dT2:dT2, R:0, P:0};
  var LMTD = (dT1 - dT2) / Math.log(dT1 / dT2);
  // R and P factors for correction factor charts
  var R_f = (Th_in - Th_out) / Math.max(Tc_out - Tc_in, 0.01);
  var P_f = (Tc_out - Tc_in) / Math.max(Th_in - Tc_in, 0.01);
  return {LMTD:LMTD, dT1:dT1, dT2:dT2, R:R_f, P:P_f};
}

var UNIT_CATS = [
  {id:"temp",cn:"\u6e29\u5ea6",en:"Temperature",units:[
    {id:"K",lb:"K",f:function(v){return v},t:function(v){return v}},
    {id:"C",lb:"\u00b0C",f:function(v){return v+273.15},t:function(v){return v-273.15}},
    {id:"F",lb:"\u00b0F",f:function(v){return (v-32)*5/9+273.15},t:function(v){return (v-273.15)*9/5+32}},
    {id:"R",lb:"\u00b0R",f:function(v){return v/1.8},t:function(v){return v*1.8}},
  ]},
  {id:"pres",cn:"\u538b\u529b",en:"Pressure",units:[
    {id:"Pa",lb:"Pa",f:function(v){return v},t:function(v){return v}},
    {id:"kPa",lb:"kPa",f:function(v){return v*1e3},t:function(v){return v/1e3}},
    {id:"MPa",lb:"MPa",f:function(v){return v*1e6},t:function(v){return v/1e6}},
    {id:"bar",lb:"bar",f:function(v){return v*1e5},t:function(v){return v/1e5}},
    {id:"atm",lb:"atm",f:function(v){return v*101325},t:function(v){return v/101325}},
    {id:"psi",lb:"psi",f:function(v){return v*6894.757},t:function(v){return v/6894.757}},
    {id:"mmHg",lb:"mmHg",f:function(v){return v*133.322},t:function(v){return v/133.322}},
    {id:"inH2O",lb:"inH2O",f:function(v){return v*249.089},t:function(v){return v/249.089}},
  ]},
  {id:"flow",cn:"\u6d41\u91cf",en:"Flow Rate",units:[
    {id:"m3h",lb:"m\u00b3/h",f:function(v){return v},t:function(v){return v}},
    {id:"m3s",lb:"m\u00b3/s",f:function(v){return v*3600},t:function(v){return v/3600}},
    {id:"Lmin",lb:"L/min",f:function(v){return v/1000*60},t:function(v){return v*1000/60}},
    {id:"gpm",lb:"gpm(US)",f:function(v){return v/3.78541*60/1000},t:function(v){return v*3.78541/60*1000}},
    {id:"bpd",lb:"bbl/day",f:function(v){return v*0.158987/24},t:function(v){return v/0.158987*24}},
  ]},
  {id:"energy",cn:"\u80fd\u91cf",en:"Energy",units:[
    {id:"J",lb:"J",f:function(v){return v},t:function(v){return v}},
    {id:"kJ",lb:"kJ",f:function(v){return v*1e3},t:function(v){return v/1e3}},
    {id:"cal",lb:"cal",f:function(v){return v*4.184},t:function(v){return v/4.184}},
    {id:"kcal",lb:"kcal",f:function(v){return v*4184},t:function(v){return v/4184}},
    {id:"BTU",lb:"BTU",f:function(v){return v*1055.06},t:function(v){return v/1055.06}},
    {id:"kWh",lb:"kWh",f:function(v){return v*3.6e6},t:function(v){return v/3.6e6}},
    {id:"MMBtu",lb:"MMBtu",f:function(v){return v*1.05506e9},t:function(v){return v/1.05506e9}},
    {id:"therm",lb:"therm",f:function(v){return v*1.05506e8},t:function(v){return v/1.05506e8}},
    {id:"toe",lb:"toe",f:function(v){return v*4.1868e10},t:function(v){return v/4.1868e10}},
  ]},
  {id:"length",cn:"\u957f\u5ea6",en:"Length",units:[
    {id:"m",lb:"m",f:function(v){return v},t:function(v){return v}},
    {id:"mm",lb:"mm",f:function(v){return v/1000},t:function(v){return v*1000}},
    {id:"cm",lb:"cm",f:function(v){return v/100},t:function(v){return v*100}},
    {id:"in",lb:"inch",f:function(v){return v*0.0254},t:function(v){return v/0.0254}},
    {id:"ft",lb:"ft",f:function(v){return v*0.3048},t:function(v){return v/0.3048}},
  ]},
  {id:"mass",cn:"\u8d28\u91cf",en:"Mass",units:[
    {id:"kg",lb:"kg",f:function(v){return v},t:function(v){return v}},
    {id:"g",lb:"g",f:function(v){return v/1000},t:function(v){return v*1000}},
    {id:"lb",lb:"lb",f:function(v){return v*0.453592},t:function(v){return v/0.453592}},
    {id:"ton",lb:"t(metric)",f:function(v){return v*1000},t:function(v){return v/1000}},
    {id:"ston",lb:"short ton",f:function(v){return v*907.185},t:function(v){return v/907.185}},
    {id:"lton",lb:"long ton",f:function(v){return v*1016.05},t:function(v){return v/1016.05}},
  ]},
  {id:"massflow",cn:"质量流量",en:"Mass Flow",units:[{id:"kgh",lb:"kg/h",f:function(v){return v},t:function(v){return v}},{id:"kgs",lb:"kg/s",f:function(v){return v*3600},t:function(v){return v/3600}},{id:"th",lb:"t/h",f:function(v){return v/1000},t:function(v){return v*1000}},{id:"lbh",lb:"lb/h",f:function(v){return v/0.45359},t:function(v){return v*0.45359}}]},
  {id:"molflow",cn:"摩尔流量",en:"Molar Flow",units:[{id:"mols",lb:"mol/s",f:function(v){return v},t:function(v){return v}},{id:"kmolh",lb:"kmol/h",f:function(v){return v*3.6},t:function(v){return v/3.6}},{id:"lbmolh",lb:"lbmol/h",f:function(v){return v*3.6/0.45359},t:function(v){return v/3.6*0.45359}}]},
  {id:"visc",cn:"\u7c98\u5ea6",en:"Viscosity",units:[
    {id:"Pas",lb:"Pa\u00b7s",f:function(v){return v},t:function(v){return v}},
    {id:"mPas",lb:"mPa\u00b7s(cP)",f:function(v){return v/1000},t:function(v){return v*1000}},
    {id:"P",lb:"Poise",f:function(v){return v/10},t:function(v){return v*10}},
  ]},
];


// PIPE WALL THICKNESS - ASME B31.3-2022
var PIPE_SCH = [
  {dn:15,od:21.3,s:[["5S",1.65],["10S",2.11],["40",2.77],["80",3.73]]},
  {dn:25,od:33.4,s:[["5S",1.65],["10S",2.77],["40",3.38],["80",4.55]]},
  {dn:40,od:48.3,s:[["5S",1.65],["10S",2.77],["40",3.68],["80",5.08]]},
  {dn:50,od:60.3,s:[["5S",1.65],["10S",2.77],["40",3.91],["80",5.54]]},
  {dn:80,od:88.9,s:[["5S",2.11],["10S",3.05],["40",5.49],["80",7.62]]},
  {dn:100,od:114.3,s:[["5S",2.11],["10S",3.05],["40",6.02],["80",8.56]]},
  {dn:150,od:168.3,s:[["5S",2.77],["10S",3.40],["40",7.11],["80",10.97]]},
  {dn:200,od:219.1,s:[["5S",2.77],["10S",3.76],["40",8.18],["80",12.70]]},
  {dn:250,od:273.1,s:[["5S",3.40],["10S",4.19],["40",9.27],["80",15.09]]},
  {dn:300,od:323.9,s:[["5S",3.96],["10S",4.57],["40",9.53],["80",17.48]]},
  {dn:400,od:406.4,s:[["5S",4.19],["10",4.78],["40",12.70],["80",21.44]]},
  {dn:500,od:508.0,s:[["5S",4.78],["10",5.54],["40",15.09],["80",26.19]]},
  {dn:600,od:609.6,s:[["5S",5.54],["10",6.35],["40",17.48],["80",30.96]]},
];
var PIPE_MAT = [
  {id:"A106B",en:"A106-B [ASME]",cn:"A106-B [ASME]",S:{20:137.9,100:137.9,200:130.3,300:114.4,400:103.4},std:"ASME"},
  {id:"20G",en:"20# [GB/T 20801]",cn:"20#\u78b3\u94a2 [GB]",S:{20:130,100:130,200:122,300:112,400:98},std:"GB"},
  {id:"TP304",en:"TP304 [ASME]",cn:"TP304 [ASME]",S:{20:137.9,100:115.1,200:98.9,300:86.9,400:79.3},std:"ASME"},
  {id:"0Cr18",en:"0Cr18Ni9 [GB]",cn:"0Cr18Ni9 [GB]",S:{20:137,100:127,200:108,300:92,400:83},std:"GB"},
  {id:"TP316L",en:"TP316L [ASME]",cn:"TP316L [ASME]",S:{20:137.9,100:112.4,200:96.5,300:84.1,400:75.8},std:"ASME"},
  {id:"00Cr17",en:"00Cr17Ni14Mo2 [GB]",cn:"00Cr17Ni14Mo2 [GB]",S:{20:133,100:120,200:100,300:88,400:78},std:"GB"},
  {id:"P12",en:"A335-P12 [ASME]",cn:"A335-P12 [ASME]",S:{20:137.9,100:137.9,200:131.0,300:124.1,400:117.2,500:82.7},std:"ASME"},
  {id:"15CrMo",en:"15CrMoR [GB]",cn:"15CrMoR [GB]",S:{20:130,100:130,200:126,300:118,400:109,500:78},std:"GB"},
  {id:"Q345R",en:"Q345R [GB]",cn:"Q345R\u5bb9\u5668\u94a2 [GB]",S:{20:189,100:185,200:174,300:162,400:143},std:"GB"},
  {id:"16MnR",en:"16MnR [GB]",cn:"16MnR [GB]",S:{20:170,100:170,200:158,300:143,400:127},std:"GB"},
];
function interpS(mat,T){var ts=[],vs=[];for(var k in mat.S){ts.push(parseFloat(k));vs.push(mat.S[k])}if(T<=ts[0])return vs[0];if(T>=ts[ts.length-1])return vs[ts.length-1];for(var i=0;i<ts.length-1;i++){if(T>=ts[i]&&T<=ts[i+1])return vs[i]+(T-ts[i])/(ts[i+1]-ts[i])*(vs[i+1]-vs[i])}return vs[0]}
function calcPipeWall(P,T,matId,dnI,Ew,Yc,Cc){var mt=PIPE_MAT.find(function(m){return m.id===matId});if(!mt)return{err:"No material"};var pp=PIPE_SCH[dnI];if(!pp)return{err:"No pipe"};var S=interpS(mt,T);var tMin=P*pp.od/(2*S*Ew+2*Yc*P)+Cc;var rec=null;for(var i=0;i<pp.s.length;i++){if(pp.s[i][1]>=tMin){rec=pp.s[i];break}}return{tMin:tMin,S:S,dn:pp.dn,od:pp.od,rec:rec,all:pp.s,std:mt.std,matName:mt.cn}}

// SAFETY VALVE - API 520 Part I
// Gas: A(mm2), W(kg/h), T(K), Z(-), MW(g/mol), gamma(-), P1(kPa_a), Pb(kPa_a), Kd(-)
// Liquid: A(mm2), Q(L/min), G(SG), P1(kPa_a), P2(kPa_a), Kd(-)
function apiC(g){return 520*Math.sqrt(g*Math.pow(2/(g+1),(g+1)/(g-1)))}
function calcSVgas(W,T,Z,MW,g,P1,Pb,Kd){
  Kd=Kd||0.975;
  var rc=Math.pow(2/(g+1),g/(g-1)); // critical pressure ratio
  var Pc_abs=P1*rc;
  var cr=Pb<=Pc_abs; // critical (choked) flow
  var Cv=apiC(g);
  var Kb=1;
  if(!cr){
    // Subcritical Kb: ratio of subcritical to critical isentropic mass flux
    // Ref: API 520 Part 1, Annex F; Darby "Chemical Engineering Fluid Mechanics" 3rd ed.
    var r=Pb/P1;
    var rc2g=Math.pow(rc,2/g), rcgm1g=Math.pow(rc,(g-1)/g);
    var den=Math.sqrt(2*g/(g-1)*rc2g*(1-rcgm1g));
    var num=Math.sqrt(2*g/(g-1)*Math.pow(r,2/g)*(1-Math.pow(r,(g-1)/g)));
    Kb=(den>0)?num/den:1;
    if(Kb>1)Kb=1; // safety clamp
  }
  var A=13160*W*Math.sqrt(T*Z)/(Cv*Kd*P1*Kb*Math.sqrt(MW));
  var orf=[["D",71],["E",126],["F",198],["G",325],["H",506],["J",830],["K",1186],["L",1841],["M",2323],["N",2800],["P",4116],["Q",6452],["R",10323],["T",16774]];
  var rec=null;for(var i=0;i<orf.length;i++){if(orf[i][1]>=A){rec=orf[i];break}}
  return{A:A,cr:cr,C:Cv,Kb:Kb,Pc:Pc_abs,rec:rec}
}
function calcSVliq(Q,G,P1,P2,Kd,Kw,Kv){Kd=Kd||0.65;Kw=Kw||1;Kv=Kv||1;var dP=P1-P2;if(dP<=0)return{err:"P1>P2"};var A=11.78*Q*Math.sqrt(G/dP)/(Kd*Kw*Kv);var orf=[["D",71],["E",126],["F",198],["G",325],["H",506],["J",830],["K",1186],["L",1841],["M",2323],["N",2800],["P",4116],["Q",6452],["R",10323],["T",16774]];var rec=null;for(var i=0;i<orf.length;i++){if(orf[i][1]>=A){rec=orf[i];break}}return{A:A,rec:rec}}

// ORIFICE PLATE - ISO 5167-2:2003 Reader-Harris/Gallagher
function calcOrifice(Dmm,dmm,dP,rho,mu,tap){var D=Dmm/1000;var d=dmm/1000;var beta=d/D;if(beta<0.1||beta>0.75)return{err:"Beta 0.1-0.75"};var Ao=Math.PI/4*d*d;var Ev=1/Math.sqrt(1-Math.pow(beta,4));var L1,L2p;if(tap==="corner"){L1=0;L2p=0}else if(tap==="D"){L1=1;L2p=0.47}else{L1=25.4/(D*1000);L2p=L1}var Cd=0.6;var Qr,v,Re;for(var it=0;it<50;it++){Qr=Cd*Ev*Ao*Math.sqrt(2*dP/rho);v=Qr/(Math.PI/4*D*D);Re=rho*v*D/mu;if(Re<5000&&it===0){return{err:"Re="+Re.toFixed(0)+" < 5000 (ISO 5167-2 minimum). Increase flow rate or pipe/orifice size."}}var Ar=Math.pow(19000*beta/Re,0.8);var M2=2*L2p/(1-beta);var Cn=0.5961+0.0261*beta*beta-0.216*Math.pow(beta,8)+0.000521*Math.pow(1e6*beta/Re,0.7)+(0.0188+0.0063*Ar)*Math.pow(beta,3.5)*Math.pow(1e6/Re,0.3)+(0.043+0.080*Math.exp(-10*L1)-0.123*Math.exp(-7*L1))*(1-0.11*Ar)*Math.pow(beta,4)/(1-Math.pow(beta,4))-0.031*(M2-0.8*Math.pow(M2,1.1))*Math.pow(beta,1.3);if(D<0.07112)Cn+=0.011*(0.75-beta)*(2.8-D/0.0254);if(Math.abs(Cn-Cd)<1e-8){Cd=Cn;break}Cd=Cd+0.5*(Cn-Cd)}Qr=Cd*Ev*Ao*Math.sqrt(2*dP/rho);v=Qr/(Math.PI/4*D*D);Re=rho*v*D/mu;return{C:Cd,E:Ev,beta:beta,Q:Qr,Qh:Qr*3600,v:v,Re:Re}}

// GAS STANDARD FLOW
function gasStdFlow(Qa,T,P,Ts,Ps){return Qa*(P/Ps)*(Ts/T)}

// MIXTURE MW
function mixMW(comps,fracs,isMass){var nc=comps.length;if(isMass){var nf=[],s=0;for(var i=0;i<nc;i++){var n=fracs[i]/comps[i].MW;nf.push(n);s+=n}for(var j=0;j<nc;j++)nf[j]/=s;var mw=0;for(var k=0;k<nc;k++)mw+=nf[k]*comps[k].MW;return{MW:mw,mol:nf,mass:fracs}}var mw2=0;for(var m=0;m<nc;m++)mw2+=fracs[m]*comps[m].MW;var mf=[];for(var p=0;p<nc;p++)mf.push(fracs[p]*comps[p].MW/mw2);return{MW:mw2,mol:fracs,mass:mf}}

// ============================================================================

// INSULATION HEAT LOSS - bare & insulated pipe/vessel

var INSUL_MATS = [
  {id:"rockwool",cn:"\u5ca9\u68c9",en:"Rock wool",k:0.040,tmax:600},
  {id:"glassfiber",cn:"\u7eff\u7483\u68c9",en:"Glass fiber",k:0.035,tmax:350},
  {id:"pufoam",cn:"\u805a\u6c28\u916f",en:"PU foam",k:0.025,tmax:120},
  {id:"casil",cn:"\u7845\u9178\u9499",en:"Calcium silicate",k:0.055,tmax:1000},
  {id:"perlite",cn:"\u73e0\u5149\u7802",en:"Perlite",k:0.050,tmax:800},
  {id:"aerogel",cn:"\u6c14\u51dd\u80f6",en:"Aerogel",k:0.015,tmax:650},
  {id:"vacuum",cn:"\u771f\u7a7a\u5939\u5957",en:"Vacuum jacket",k:0.001,tmax:200},
  {id:"custom",cn:"\u81ea\u5b9a\u4e49",en:"Custom",k:0,tmax:9999},
];
function calcInsulation(D_mm, Tproc_C, Tamb_C, insThick_mm, insK, windV, emissivity) {
  var D = D_mm / 1000;
  var Tp = Tproc_C + 273.15;
  var Ta = Tamb_C + 273.15;
  var sigma = 5.67e-8;
  emissivity = emissivity || 0.9;
  windV = windV || 3;
  insK = insK || 0.04;
  var hConv = 5.7 + 3.8 * windV;
  var hRad = emissivity * sigma * (Tp * Tp + Ta * Ta) * (Tp + Ta);
  var Uo_bare = hConv + hRad;
  var Qbare = Uo_bare * Math.PI * D * (Tp - Ta);
  var Tsurf_bare = Tp;
  if (insThick_mm <= 0) return {Qbare: Qbare, Qins: Qbare, Tsurf_bare: Tsurf_bare - 273.15, Tsurf_ins: Tsurf_bare - 273.15, saving: 0, Uo_bare: Uo_bare};
  var Dins = D + 2 * insThick_mm / 1000;
  var Rins = Math.log(Dins / D) / (2 * Math.PI * insK);
  var hConvIns = 5.7 + 3.8 * windV;
  var Ro = 1 / (hConvIns * Math.PI * Dins);
  var Rtot = Rins + Ro;
  var Qins = (Tp - Ta) / Rtot;
  var Tsurf_ins = Ta + Qins * Ro;
  var saving = (1 - Qins / Qbare) * 100;
  return {Qbare: Qbare, Qins: Qins, Tsurf_bare: Tsurf_bare - 273.15, Tsurf_ins: Tsurf_ins - 273.15, saving: saving, Uo_bare: Uo_bare, Dins: Dins * 1000};
}

// WATER PHASE BOUNDARY DATA (simplified for P-T diagram)
var WATER_PHASE = {
  triple: {T: 0.01, P: 0.000612},
  critical: {T: 373.946, P: 22.064},
  getData: function() {
    var pts = [];
    // Sublimation line (solid-vapor): A -> triple point
    var sT = [-40,-35,-30,-25,-20,-15,-10,-5,0,0.01];
    var sP = [0.0000129,0.0000223,0.0000381,0.0000632,0.000103,0.000165,0.000260,0.000401,0.000611,0.000612];
    for (var i = 0; i < sT.length; i++) pts.push({T:sT[i],subl:sP[i]});
    // Vaporization line (liquid-vapor): triple point -> critical point
    var vT = [0.01,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240,260,280,300,320,340,360,373.946];
    var vP = [0.000612,0.001228,0.002339,0.004247,0.007384,0.01235,0.01995,0.03120,0.04740,0.07018,0.10142,0.19867,0.3614,0.6182,1.0028,1.5549,2.320,3.348,4.692,6.417,8.588,11.284,14.601,18.666,22.064];
    for (var j = 0; j < vT.length; j++) pts.push({T:vT[j],vap:vP[j]});
    // Melting line (solid-liquid): triple point -> upward (water: negative slope!)
    var mP = [0.000612,0.101,1,10,50,100,200,400,600];
    var mT = [0.01,0.00,-0.0075,-0.075,-0.375,-0.75,-1.5,-3.0,-4.5];
    for (var k = 0; k < mP.length; k++) pts.push({T:mT[k],melt:mP[k]});
    return pts;
  }
};
var CO2_PHASE = {
  triple: {T: -56.56, P: 0.5175},
  critical: {T: 30.98, P: 7.377},
  getData: function() {
    var pts = [];
    // Sublimation line
    var sT = [-120,-110,-100,-95,-90,-85,-80,-78.5,-75,-70,-65,-60,-56.56];
    var sP = [0.00028,0.00144,0.00672,0.01342,0.02533,0.04540,0.07760,0.10133,0.14,0.21,0.327,0.464,0.5175];
    for (var i = 0; i < sT.length; i++) pts.push({T:sT[i],subl:sP[i]});
    // Vaporization line
    var vT = [-56.56,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,30.98];
    var vP = [0.5175,0.683,0.852,1.005,1.213,1.464,1.753,2.085,2.466,2.903,3.203,3.485,3.970,4.502,5.086,5.729,6.434,7.11,7.377];
    for (var j = 0; j < vT.length; j++) pts.push({T:vT[j],vap:vP[j]});
    // Melting line (positive slope for CO2)
    var mP = [0.5175,5,10,20,50,100,200];
    var mT = [-56.56,-56.0,-55.4,-54.2,-51.5,-46.5,-36.5];
    for (var k = 0; k < mP.length; k++) pts.push({T:mT[k],melt:mP[k]});
    return pts;
  }
};

// 4-20mA SIGNAL CONVERSION
function maConvert(pv, lo, hi, maLo, maHi) {
  maLo = maLo || 4; maHi = maHi || 20;
  var pct = (pv - lo) / (hi - lo);
  var ma = maLo + pct * (maHi - maLo);
  return {pct: pct * 100, ma: ma, pv: pv, lo: lo, hi: hi};
}
function maReverse(ma, lo, hi, maLo, maHi) {
  maLo = maLo || 4; maHi = maHi || 20;
  var pct = (ma - maLo) / (maHi - maLo);
  var pv = lo + pct * (hi - lo);
  return {pct: pct * 100, ma: ma, pv: pv, lo: lo, hi: hi};
}

// PIPE VELOCITY RECOMMENDATIONS (m/s)
var PIPE_VEL = [
  {fluid:"Water (suction)",cn:"\u6c34(\u5438\u5165)",lo:0.5,hi:1.5,typ:1.0},
  {fluid:"Water (discharge)",cn:"\u6c34(\u6392\u51fa)",lo:1.5,hi:3.0,typ:2.0},
  {fluid:"Steam (low P <0.5MPa)",cn:"\u84b8\u6c7d(\u4f4eP)",lo:15,hi:30,typ:20},
  {fluid:"Steam (high P >1MPa)",cn:"\u84b8\u6c7d(\u9ad8P)",lo:25,hi:50,typ:35},
  {fluid:"Air (compressed)",cn:"\u538b\u7f29\u7a7a\u6c14",lo:8,hi:20,typ:15},
  {fluid:"Natural gas",cn:"\u5929\u7136\u6c14",lo:10,hi:25,typ:15},
  {fluid:"H2 gas",cn:"\u6c22\u6c14",lo:15,hi:40,typ:25},
  {fluid:"N2 gas",cn:"\u6c2e\u6c14",lo:10,hi:25,typ:18},
  {fluid:"Organic liquid",cn:"\u6709\u673a\u6db2\u4f53",lo:0.5,hi:2.0,typ:1.2},
  {fluid:"Oil (viscous)",cn:"\u91cd\u6cb9(\u9ad8\u7c98)",lo:0.3,hi:1.0,typ:0.5},
  {fluid:"Refrigerant (liquid)",cn:"\u5236\u51b7\u5242(\u6db2)",lo:0.5,hi:1.5,typ:1.0},
  {fluid:"Refrigerant (vapor)",cn:"\u5236\u51b7\u5242(\u6c14)",lo:5,hi:15,typ:10},
  {fluid:"Gravity drain",cn:"\u91cd\u529b\u6392\u6c34",lo:0.3,hi:0.6,typ:0.5},
  {fluid:"Flue gas",cn:"\u70df\u6c14",lo:10,hi:20,typ:15},
];

// HEAT EXCHANGER AREA ESTIMATION Q=UA*LMTD
var HX_U = [
  {hot:"Water",cold:"Water",cn:"\u6c34-\u6c34",lo:800,hi:2500,typ:1500},
  {hot:"Steam",cold:"Water",cn:"\u84b8\u6c7d-\u6c34",lo:1000,hi:4000,typ:2500},
  {hot:"Steam",cold:"Organic",cn:"\u84b8\u6c7d-\u6709\u673a",lo:300,hi:1200,typ:700},
  {hot:"Steam",cold:"Gas",cn:"\u84b8\u6c7d-\u6c14\u4f53",lo:20,hi:300,typ:100},
  {hot:"Hot oil",cold:"Water",cn:"\u5bfc\u70ed\u6cb9-\u6c34",lo:300,hi:1000,typ:600},
  {hot:"Hot oil",cold:"Organic",cn:"\u5bfc\u70ed\u6cb9-\u6709\u673a",lo:100,hi:500,typ:300},
  {hot:"Gas",cold:"Gas",cn:"\u6c14-\u6c14",lo:10,hi:50,typ:25},
  {hot:"Gas",cold:"Water",cn:"\u6c14\u4f53-\u6c34",lo:20,hi:300,typ:100},
  {hot:"Condensing steam",cold:"Boiling water",cn:"\u51dd\u6c7d-\u6cb8\u6c34",lo:2000,hi:5000,typ:3500},
  {hot:"Organic cond.",cold:"Water",cn:"\u6709\u673a\u51dd-\u6c34",lo:300,hi:1200,typ:700},
];
function calcHXarea(Q_kW, U, LMTD_K) {
  if (U <= 0 || LMTD_K <= 0) return {err: "U,LMTD>0"};
  var A = Q_kW * 1000 / (U * LMTD_K);
  return {A: A, Q: Q_kW, U: U, LMTD: LMTD_K};
}

// THERMOCOUPLE mV->T LOOKUP (ITS-90 simplified polynomial)
// Thermocouple types: K, B, S (ITS-90) + Pt100, Pt1000 (RTD)
// B-type: PtRh30%-PtRh6%, range 0-1820°C, ITS-90 coefficients
// S-type: PtRh10%-Pt, range -50-1768°C, ITS-90 coefficients
// K-type: NiCr-NiAl, range -200-1372°C
var TC_TYPES = {
  Pt100: {name:"Pt100 (RTD)",range:[-200,850],
    toOhm:function(T){if(T>=0)return 100*(1+3.9083e-3*T-5.775e-7*T*T);return 100*(1+3.9083e-3*T-5.775e-7*T*T-4.183e-12*Math.pow(T-100,3)*T)},
    toT:function(R){var a=3.9083e-3,b=-5.775e-7;return(-a+Math.sqrt(a*a-4*b*(1-R/100)))/(2*b)}
  },
  Pt1000: {name:"Pt1000 (RTD)",range:[-200,850],
    toOhm:function(T){if(T>=0)return 1000*(1+3.9083e-3*T-5.775e-7*T*T);return 1000*(1+3.9083e-3*T-5.775e-7*T*T-4.183e-12*Math.pow(T-100,3)*T)},
    toT:function(R){var a=3.9083e-3,b=-5.775e-7;return(-a+Math.sqrt(a*a-4*b*(1-R/1000)))/(2*b)}
  },
  K: {name:"K 型 (NiCr-NiAl)",range:[-200,1372],
    toMV:function(T){if(T<0)return 0.03945+T*(0.03969+T*(4.21e-5+T*2.17e-8));return-0.0176+T*(0.03974+T*(2.95e-5+T*(-3.33e-8+T*1.24e-11)))},
    toT:function(mV){if(mV<0)return mV*(25.08+mV*(0.072+mV*0.0022));return mV*(25.08+mV*(-0.608+mV*(0.0145+mV*(-0.000178))))}
  },
  B: {name:"B 型 (PtRh30-PtRh6)",range:[0,1820],
    // Polynomial fitted to NIST ITS-90 B-type tabulated data. Max error: <0.002 mV
    toMV:function(T){
      if(T<0||T>1820)return NaN;
      if(T<=630){
        var c=[4.32900433e-6,-2.45227273e-4,5.82007576e-6,-7.04545455e-10,-3.78787879e-14,0];
        var E=0;for(var i=0;i<c.length;i++)E+=c[i]*Math.pow(T,i);return E;
      } else {
        var d=[-6.23505199e-1,2.31957556e-3,1.91160520e-6,1.93331486e-9,-7.07683194e-13];
        var E=0;for(var i=0;i<d.length;i++)E+=d[i]*Math.pow(T,i);return E;
      }
    },
    toT:function(mV){
      if(mV<0||mV>13.82)return NaN;
      var T=mV>2?900:300;
      for(var it=0;it<40;it++){
        var f=TC_TYPES.B.toMV(T)-mV;
        var df=(TC_TYPES.B.toMV(T+0.5)-TC_TYPES.B.toMV(T-0.5))/1.0;
        if(!isFinite(df)||Math.abs(df)<1e-10)break;
        var step=f/df;if(Math.abs(step)>100)step=100*(step>0?1:-1);
        T-=step;T=Math.max(0,Math.min(1820,T));
        if(Math.abs(step)<0.01)break;
      }
      return T;
    }
  },
  S: {name:"S 型 (PtRh10-Pt)",range:[0,1768],
    // Polynomial fitted to NIST ITS-90 S-type tabulated data. Max error: <0.006 mV
    toMV:function(T){
      if(T<0||T>1768)return NaN;
      if(T<=1064){
        var c=[-4.14335664e-4,5.51046300e-3,1.11215691e-5,-1.56412442e-8,1.23302739e-11,-3.73397436e-15];
        var E=0;for(var i=0;i<c.length;i++)E+=c[i]*Math.pow(T,i);return E;
      } else {
        var d=[3.18991670e0,-7.96086287e-4,9.55869097e-6,-2.34560382e-9];
        var E=0;for(var i=0;i<d.length;i++)E+=d[i]*Math.pow(T,i);return E;
      }
    },
    toT:function(mV){
      if(mV<0||mV>18.693)return NaN;
      // Better initial guess: linear interpolation between key points
      var T = mV<9.587 ? mV/9.587*1000 : 1000 + (mV-9.587)/(18.693-9.587)*768;
      T = Math.max(0, Math.min(1768, T));
      for(var it=0;it<50;it++){
        var f=TC_TYPES.S.toMV(T)-mV;
        var df=(TC_TYPES.S.toMV(Math.min(T+0.5,1768))-TC_TYPES.S.toMV(Math.max(T-0.5,0)))/1.0;
        if(!isFinite(df)||Math.abs(df)<1e-10)break;
        var step=f/df;if(Math.abs(step)>100)step=100*(step>0?1:-1);
        T-=step;T=Math.max(0,Math.min(1768,T));
        if(Math.abs(step)<0.005)break;
      }
      return T;
    }
  },
};

;

// NPSH CHECK
function calcNPSH(Ps_kPa, Pv_kPa, hs_m, hf_m, rho) {
  rho = rho || 1000;
  var g = 9.81;
  var NPSHa = (Ps_kPa - Pv_kPa) * 1000 / (rho * g) + hs_m - hf_m;
  return {NPSHa: NPSHa, Ps: Ps_kPa, Pv: Pv_kPa, hs: hs_m, hf: hf_m};
}



// MULTI-STAGE COMPRESSOR
function calcMultiComp(P1,P2,T1,eta,gam,mw,flow,nst,Tcool_C){
  var TcoolK = Tcool_C + 273.15;
  if(nst<=1){
    var T2s=T1*Math.pow(P2/P1,(gam-1)/gam);
    var T2=T1+(T2s-T1)/eta;
    var Ws=flow/3600*(gam/(gam-1))*(R_GAS/mw*1000)*T1*(Math.pow(P2/P1,(gam-1)/gam)-1)/eta;
    return {T2s:T2s,T2:T2,W:Ws/1000,nst:1,rps:P2/P1,stages:[{pin:P1/1e6,pout:P2/1e6,tin:T1-273.15,t2s:T2s-273.15,t2:T2-273.15,w:Ws/1000}]};
  }
  var rps=Math.pow(P2/P1,1/nst);
  var stages=[],totalW=0,Tin=T1;
  for(var si=0;si<nst;si++){
    var Pin=P1*Math.pow(rps,si),Pout=Pin*rps;
    var T2si=Tin*Math.pow(rps,(gam-1)/gam);
    var T2i=Tin+(T2si-Tin)/eta;
    var Wi=flow/3600*(gam/(gam-1))*(R_GAS/mw*1000)*Tin*(Math.pow(rps,(gam-1)/gam)-1)/eta;
    stages.push({pin:Pin/1e6,pout:Pout/1e6,tin:Tin-273.15,t2s:T2si-273.15,t2:T2i-273.15,w:Wi/1000});
    totalW+=Wi;
    Tin=si<nst-1?TcoolK:T2i; // Use TcoolK (in Kelvin) for intercooler
  }
  return {nst:nst,stages:stages,W:totalW/1000,rps:rps,T2:stages[nst-1].t2+273.15,T2s:stages[nst-1].t2s+273.15};
}

// HUMIDITY CALCULATOR - Sonntag (1990) saturation pressure correlations
// Valid range: -150°C to +100°C (ice curve extrapolated below -100°C)
// Sonntag (1990) saturation pressure over liquid water [Pa→kPa]
// Coefficients verified: at 25°C gives 3.169 kPa, at 100°C gives 101.3 kPa
function satPressWater(T_C) {
  var T = T_C + 273.15;
  var lnP = -6096.9385/T + 21.2409642 - 0.02711193*T + 1.673952e-5*T*T + 2.433502*Math.log(T);
  return Math.exp(lnP) / 1000; // Pa → kPa
}
// Sonntag (1990) saturation pressure over ice [Pa→kPa]
// Coefficients verified: at 0°C gives 0.611 kPa, at -40°C gives 0.01284 kPa
function satPressIce(T_C) {
  var T = T_C + 273.15;
  var lnP = -6024.5282/T + 29.32707 + 0.010613868*T - 1.3198825e-5*T*T - 0.49382577*Math.log(T);
  return Math.exp(lnP) / 1000; // Pa → kPa
}
function satPressBest(T_C) {
  return T_C <= 0 ? satPressIce(T_C) : satPressWater(T_C);
}
function dewPointFromPw(Pw_kPa) {
  if (Pw_kPa <= 0) return -273.15;
  var a=17.625, b=243.04;
  var Td = b*Math.log(Pw_kPa/0.61078)/(a-Math.log(Pw_kPa/0.61078));
  for (var iter=0; iter<50; iter++) {
    var Ps = satPressBest(Td);
    var diff = Ps - Pw_kPa;
    if (Math.abs(diff) < 1e-8) break;
    var Ps2 = satPressBest(Td + 0.01);
    var dPdT = (Ps2 - Ps) / 0.01;
    if (Math.abs(dPdT) < 1e-15) break;
    var step = -diff / dPdT;
    if (Math.abs(step) > 10) step = 10*(step>0?1:-1);
    Td += step;
    if (Td < -150) Td = -150;
    if (Td > 100) Td = 100;
  }
  return Td;
}
function frostPointFromPw(Pw_kPa) {
  if (Pw_kPa <= 0) return -273.15;
  var Tf = -40;
  for (var iter=0; iter<50; iter++) {
    var Ps = satPressIce(Tf);
    var diff = Ps - Pw_kPa;
    if (Math.abs(diff) < 1e-8) break;
    var Ps2 = satPressIce(Tf + 0.01);
    var dPdT = (Ps2 - Ps) / 0.01;
    if (Math.abs(dPdT) < 1e-15) break;
    var step = -diff / dPdT;
    if (Math.abs(step) > 10) step = 10*(step>0?1:-1);
    Tf += step;
    if (Tf < -150) Tf = -150;
    if (Tf > 0) Tf = 0;
  }
  return Tf;
}
function calcHumidity(Tdb, RH, Patm) {
  var Psat = satPressBest(Tdb);
  var Pw = Psat * RH / 100;
  var Tdew = dewPointFromPw(Pw);
  var Tfrost = Tdew < 0 ? frostPointFromPw(Pw) : null;
  var W = 0.62198 * Pw / (Patm - Pw);
  var Twet = Tdb*Math.atan(0.151977*Math.pow(RH+8.3137,0.5))+Math.atan(Tdb+RH)-Math.atan(RH-1.6763)+0.00391838*Math.pow(RH,1.5)*Math.atan(0.023101*RH)-4.686035;
  return {Psat:Psat, Pw:Pw, Tdew:Tdew, Tfrost:Tfrost, W:W*1000, Twet:Twet, ppm:Pw/Patm*1e6};
}

// B36.10M PIPE SCHEDULE
var B3610=[{nps:"1/2",dn:15,od:21.3,sch:{STD:2.77,XS:3.73,"160":4.78,XXS:7.47}},{nps:"3/4",dn:20,od:26.7,sch:{STD:2.87,XS:3.91,"160":5.56,XXS:7.82}},{nps:"1",dn:25,od:33.4,sch:{STD:3.38,XS:4.55,"160":6.35,XXS:9.09}},{nps:"1-1/2",dn:40,od:48.3,sch:{STD:3.68,XS:5.08,"160":7.14,XXS:10.15}},{nps:"2",dn:50,od:60.3,sch:{STD:3.91,XS:5.54,"160":8.74,XXS:11.07}},{nps:"3",dn:80,od:88.9,sch:{STD:5.49,XS:7.62,"160":11.13,XXS:15.24}},{nps:"4",dn:100,od:114.3,sch:{STD:6.02,XS:8.56,"160":13.49,XXS:17.12}},{nps:"6",dn:150,od:168.3,sch:{STD:7.11,XS:10.97,"160":18.26,XXS:21.95}},{nps:"8",dn:200,od:219.1,sch:{STD:8.18,XS:12.70,"160":23.01,XXS:22.23}},{nps:"10",dn:250,od:273.1,sch:{STD:9.27,XS:12.70,"160":28.58,XXS:25.40}},{nps:"12",dn:300,od:323.8,sch:{STD:9.53,XS:12.70,"160":33.32,XXS:25.40}},{nps:"14",dn:350,od:355.6,sch:{STD:9.53,XS:12.70,"160":35.71}},{nps:"16",dn:400,od:406.4,sch:{STD:9.53,XS:12.70,"160":40.49}},{nps:"20",dn:500,od:508.0,sch:{STD:9.53,XS:12.70,"160":50.01}},{nps:"24",dn:600,od:609.6,sch:{STD:9.53,XS:12.70,"160":59.54}}];

// ============================================================================
// MIXTURE PROPERTIES: Viscosity + Density
// ============================================================================
// Gas viscosity: Herning-Zipperer mixing rule
// mu_mix = sum(yi*mui*sqrt(MWi)) / sum(yi*sqrt(MWi))
function gasViscosityLucas(T, comp) {
  // Lucas method for pure component gas viscosity [Pa·s]
  // Valid for reduced temperature 0.2 < Tr < 100
  var Tc=comp.Tc, Pc=comp.Pc, MW=comp.MW, w=comp.w;
  var Tr=T/Tc;
  // Lucas correlating parameter Xi
  var xi = 0.176*Math.pow(Tc/Math.pow(MW,3)/Math.pow(Pc/1e5,4), 1/6); // [microP^-1]
  var FP0, FQ0;
  // Polarity factor FP0
  var mu_r = 52.46*w*Pc/1e5/(Tc*Tc); // reduced dipole moment approximation
  if (mu_r < 0.022) FP0=1;
  else if (mu_r < 0.075) FP0=1+30.55*Math.pow(0.292-w,1.72);
  else FP0=1+30.55*Math.pow(0.292-w,1.72)*Math.abs(0.96+0.1*(Tr-0.7));
  FQ0=1; // quantum correction (1 for most components except H2, He)
  if (comp.id==="H2") FQ0=1.22*Math.pow(1/Tr,0.15)*(1+0.00385*Math.pow((Tr-12)*(Tr-12)/MW,1/8)*Math.sign(Tr-12));
  if (comp.id==="He") FQ0=1.22*Math.pow(1/Tr,0.15);
  // Low-pressure viscosity
  var Z1 = (0.807*Math.pow(Tr,0.618) - 0.357*Math.exp(-0.449*Tr)
    + 0.340*Math.exp(-4.058*Tr) + 0.018) * FP0 * FQ0;
  return Z1 / xi * 1e-7; // Pa·s
}

function gasMixViscosity(T, y, clist) {
  // Herning-Zipperer: mu_mix = sum(yi*mui*sqrt(MWi)) / sum(yi*sqrt(MWi))
  var num=0, den=0;
  for (var i=0;i<clist.length;i++) {
    var mu=gasViscosityLucas(T, clist[i]);
    var sq=Math.sqrt(clist[i].MW);
    num+=y[i]*mu*sq; den+=y[i]*sq;
  }
  return den>0 ? num/den : 0;
}

function liqViscosityAndrade(T, comp) {
  // Andrade equation: ln(mu) = A + B/T
  // Calibrated at Tb (mu_b ~ 1 mPa·s estimate) and Tc (mu_c ~ 0.2 mPa·s)
  // More reliable: use Letsou-Stiel at Tr<1
  var Tc=comp.Tc, Tr=T/Tc, w=comp.w;
  if (Tr>=1) return null; // supercritical - no liquid
  // Letsou-Stiel: xi*mu = [xi*mu]^(0) + w*[xi*mu]^(1)
  var xi = Math.pow(Tc, 1/6) / (Math.pow(comp.MW, 0.5) * Math.pow(comp.Pc/1e5, 2/3));
  var xmu0 = (1.5174 - 2.135*Tr + 0.75*Tr*Tr) * 1e-5;
  var xmu1 = (4.2552 - 7.674*Tr + 3.4*Tr*Tr) * 1e-5;
  var mu = Math.max((xmu0 + w*xmu1)/xi, 1e-7);
  return mu; // Pa·s
}

function liqMixViscosity(T, x, clist) {
  // Check for water + hydrocarbon system
  var hasWater=false, hasHC=false;
  for (var i=0;i<clist.length;i++) {
    if (clist[i].id==="H2O") hasWater=true;
    if (["CH4","C2H6","C2H4","C3H8","C3H6","nC4","iC4","Benzene","Toluene","pXylene","CyC6"].indexOf(clist[i].id)>=0) hasHC=true;
  }
  var muArr=[], warnLLE=hasWater&&hasHC;
  for (var i=0;i<clist.length;i++) {
    var mu=liqViscosityAndrade(T, clist[i]);
    muArr.push(mu||1e-3);
  }
  // Lobe: ln(mu_mix) = sum(xi*ln(mui))
  var lnMu=0;
  for (var i=0;i<clist.length;i++) lnMu+=x[i]*Math.log(muArr[i]);
  return {mu: Math.exp(lnMu), warn: warnLLE};
}

function liqMixDensity(T, x, clist) {
  // Ideal volume mixing with DIPPR 105 or Rackett for each component
  var vMix=0, MWmix=0;
  for (var i=0;i<clist.length;i++) {
    var c=clist[i], Tc=c.Tc, Tr=T/Tc;
    if (Tr>=1) return null; // above critical - no liquid density
    var rho_i; // mol/m3
    var ld=LIQDEN[c.id];
    if (ld) {
      // DIPPR 105: rho[mol/m3] = A/B^(1+(1-Tr)^D)
      rho_i = ld.A / Math.pow(ld.B, 1+Math.pow(1-Tr, ld.D));
      rho_i *= 1000; // kmol/m3 -> mol/m3
    } else {
      // Rackett: V = Vc * Zc^((1-Tr)^(2/7))
      // Approximate Zc = Pc*Vc/(R*Tc)
      if (!c.Vc) return null;
      var Zc = c.Pc * c.Vc / (R_GAS * c.Tc);
      var Vm = c.Vc * Math.pow(Zc, Math.pow(1-Tr, 2/7)); // m3/mol
      rho_i = 1/Vm; // mol/m3
    }
    var Vm_i = 1/rho_i; // m3/mol
    vMix += x[i] * Vm_i;
    MWmix += x[i] * c.MW;
  }
  if (vMix<=0) return null;
  return (MWmix/1000) / vMix; // kg/m3
}

function gasMixDensity(T, P, y, clist) {
  // PR EOS Z factor for gas phase
  var mx=mixPR(T, y, clist);
  var A=mx.am*P/(R_GAS*R_GAS*T*T), B=mx.bm*P/(R_GAS*T);
  var roots=prRoots(A,B);
  if(!roots||roots.length===0) return null;
  var Z=roots[roots.length-1]; // vapor root
  var MWmix=0; for(var i=0;i<clist.length;i++) MWmix+=y[i]*clist[i].MW;
  return P*MWmix/(1000*Z*R_GAS*T); // kg/m3
}

// ============================================================================
// BURSTING DISC (爆破片) CALCULATION
// GB/T 567.2-2012 and ISO 4126-2:2003
// ============================================================================
// Gas/vapor: A = W/(C*Kd*P0) * sqrt(TZ/M)  [GB/T 567.2]
//            A = W/(0.9*Kd*P0*C) * sqrt(TZ/M)  [ISO 4126-2]
// Liquid:    A = Q/(Kd*Kw*Kv) * sqrt(G/dP)
// A in mm2, W in kg/h, P0 in kPa(a), T in K, M in g/mol
// C = gas flow constant = 520*sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)))

function calcBurstDiscGas(W, T, Z, MW, gamma, P0, std, Kd) {
  // W[kg/h], T[K], Z[-], MW[g/mol], P0[kPa_abs], std:"GB"or"ISO", Kd[-]
  Kd = Kd || 0.62;
  var C = 520 * Math.sqrt(gamma * Math.pow(2/(gamma+1), (gamma+1)/(gamma-1)));
  var A_mm2;
  if (std === "ISO") {
    // ISO 4126-2: A = W/(0.9*Kd*P0*C) * sqrt(T*Z/M) * 1e6/3600 unit factor
    // Full SI derivation gives: A[mm2] = (W[kg/h]*1000/3600) / (0.9*Kd*P0[Pa]*C/sqrt(R)) * sqrt(T*Z*R/M)
    // Simplified with empirical constant consistent with API 520 SI:
    A_mm2 = 13160 * W * Math.sqrt(T*Z) / (0.9 * Kd * C * P0 * Math.sqrt(MW));
  } else {
    // GB/T 567.2-2012 Eq: A[mm2] = 13160*W*sqrt(TZ) / (C*Kd*P0*sqrt(M))
    A_mm2 = 13160 * W * Math.sqrt(T*Z) / (C * Kd * P0 * Math.sqrt(MW));
  }
  return {A: A_mm2, C: C, Kd: Kd, std: std};
}

function calcBurstDiscLiq(Q, G, P1_kPa, P2_kPa, std, Kd, Kw, Kv) {
  // Q[L/min], G[-], P1[kPa_abs], P2[kPa_abs]
  Kd = Kd || 0.62; Kw = Kw || 1.0; Kv = Kv || 1.0;
  var dP = P1_kPa - P2_kPa;
  if (dP <= 0) return {err: "P1 必须大于 P2"};
  var factor = std === "ISO" ? 11.78 : 11.78; // same numeric constant for both
  var A_mm2 = factor * Q * Math.sqrt(G/dP) / (Kd * Kw * Kv);
  return {A: A_mm2, Kd: Kd, Kw: Kw, Kv: Kv, dP: dP, std: std};
}

// UI LAYER - Modern Clean Industrial Design
// ============================================================================


// CONTROL VALVE Cv - IEC 60534 / ISA S75
// Liquid: Cv = Q*sqrt(G/(P1-P2)) / N1
// Gas: Cv = W/(N8*Fp*Y*sqrt(x*M*P1*rho1))
var CV_N1 = 0.0865; // m3/h, kPa
var CV_N8 = 34.7;   // IEC 60534-2-1: W[kg/h],P[kPa],T[K],M[g/mol] (N8 for sqrt(x*M*P1/(T*Z)) form)
function calcCvLiq(Q_m3h, G, P1_kPa, P2_kPa) {
  // IEC 60534-2-1 liquid Cv (no choked correction - per engineering practice)
  // Cv = Q[m3/h] / N1 * sqrt(G / dP[kPa])   N1=0.0865
  var dP = P1_kPa - P2_kPa;
  if (dP <= 0) return {err: "P1必须大于P2 / P1 must > P2"};
  var Cv = Q_m3h / CV_N1 * Math.sqrt(G / dP);
  return {Cv: Cv, dP: dP};
}
function calcCvGas(W_kgh, P1_kPa, P2_kPa, T_K, MW, gamma, Z, xT) {
  // IEC 60534-2-1 gas/vapor Cv (mass flow basis)
  // Cv = W / (N8 * Y * sqrt(xeff * M * P1 / (T * Z)))   N8=34.7
  Z = Z || 1; gamma = gamma || 1.4; xT = xT || 0.7;
  var x = (P1_kPa - P2_kPa) / P1_kPa;
  var Fk = gamma / 1.4;
  var xeff = Math.min(x, Fk * xT);
  var choked = x >= Fk * xT;
  var Y = 1 - xeff / (3 * Fk * xT);
  var Cv = W_kgh / (CV_N8 * Y * Math.sqrt(xeff * MW * P1_kPa / (T_K * Z)));
  return {Cv: Cv, x: x, xeff: xeff, Y: Y, Fk: Fk, choked: choked, regime: choked ? "阻塞/Choked" : "正常/Normal"};
}

// HEATING VALUE - HHV/LHV for common gas components
// Units: kJ/mol (divide by MW for kJ/g, multiply density for kJ/Nm3)
var HV_DATA = {
  "H2":  {MW:2.016, HHV:285.8, LHV:241.8},
  "CO":  {MW:28.01, HHV:283.0, LHV:283.0},
  "CH4": {MW:16.04, HHV:890.3, LHV:802.3},
  "C2H6":{MW:30.07, HHV:1560.7,LHV:1428.6},
  "C2H4":{MW:28.05, HHV:1411.2,LHV:1323.1},
  "C3H8":{MW:44.10, HHV:2220.0,LHV:2043.1},
  "C3H6":{MW:42.08, HHV:2058.0,LHV:1926.0},
  "nC4": {MW:58.12, HHV:2878.5,LHV:2657.3},
  "iC4": {MW:58.12, HHV:2869.0,LHV:2648.0},
  "H2S": {MW:34.08, HHV:562.0, LHV:518.0},
  "NH3": {MW:17.03, HHV:382.6, LHV:316.8},
  "MeOH":{MW:32.04, HHV:726.1, LHV:638.1},
  "EtOH":{MW:46.07, HHV:1366.8,LHV:1234.8},
};
function calcHeatVal(compIds, fracs) {
  var hhv = 0, lhv = 0, mw = 0;
  for (var i = 0; i < compIds.length; i++) {
    var d = HV_DATA[compIds[i]];
    if (!d) continue;
    hhv += fracs[i] * d.HHV;
    lhv += fracs[i] * d.LHV;
    mw += fracs[i] * d.MW;
  }
  if (mw <= 0) return {err: "No combustible components"};
  var rhoSTP = mw / 22.414; // kg/Nm3
  var hhv_kg = hhv * 1000 / mw; // kJ/kg
  var lhv_kg = lhv * 1000 / mw;
  var hhv_nm3 = hhv * 1000 / 22.414; // kJ/Nm3
  var lhv_nm3 = lhv * 1000 / 22.414;
  var rhoAir = 28.97 / 22.414;
  var wobbe = hhv_nm3 / Math.sqrt(rhoSTP / rhoAir);
  return {hhv_mol: hhv, lhv_mol: lhv, hhv_kg: hhv_kg, lhv_kg: lhv_kg, hhv_nm3: hhv_nm3, lhv_nm3: lhv_nm3, wobbe: wobbe, MW: mw, rho: rhoSTP};
}

// VESSEL VOLUME & LEVEL - horizontal/vertical cylinder + heads
// Horizontal cylinder partial volume
function horizCylVol(D, L, h) {
  if (h <= 0) return 0; if (h >= D) return Math.PI / 4 * D * D * L;
  var R = D / 2;
  var A = R * R * Math.acos((R - h) / R) - (R - h) * Math.sqrt(2 * R * h - h * h);
  return A * L;
}
// Ellipsoidal head (2:1) partial volume
function ellipHead(D, h) {
  var R = D / 2; var a = R / 2; // semi-axis for 2:1 elliptical
  if (h <= 0) return 0; if (h >= D) return 2 / 3 * Math.PI * R * R * a;
  return Math.PI * a / (R * R) * (R * R * h * h / 2 - h * h * h * h / (4 * R * R) * R * R);
}
function calcVessel(D_m, L_m, h_m, headType, orient) {
  if (D_m <= 0 || L_m <= 0) return {err: "D,L > 0"};
  var hh = Math.max(0, Math.min(h_m, D_m));
  var Vcyl, Vhead, Vtotal, VcylFull, VheadFull, VtotalFull;
  if (orient === "horiz") {
    Vcyl = horizCylVol(D_m, L_m, hh);
    VcylFull = Math.PI / 4 * D_m * D_m * L_m;
    if (headType === "flat") { Vhead = 0; VheadFull = 0; }
    else {
      // 2:1 elliptical head approximate: for horizontal, use same cross-section
      var headLen = D_m / 4; // depth of 2:1 elliptical head
      Vhead = 2 * horizCylVol(D_m, headLen, hh) * 0.667;
      VheadFull = 2 * Math.PI / 4 * D_m * D_m * headLen * 0.667;
    }
  } else {
    // Vertical: simple cylinder
    var hCyl = Math.max(0, Math.min(hh, L_m));
    Vcyl = Math.PI / 4 * D_m * D_m * hCyl;
    VcylFull = Math.PI / 4 * D_m * D_m * L_m;
    if (headType === "flat") { Vhead = 0; VheadFull = 0; }
    else {
      var headH = D_m / 4;
      VheadFull = 2 / 3 * Math.PI * Math.pow(D_m / 2, 2) * headH;
      // Bottom head fills first
      var hBot = Math.max(0, Math.min(hh, headH));
      var hTop = Math.max(0, hh - L_m - headH);
      Vhead = VheadFull * Math.pow(hBot / headH, 2) * (3 - 2 * hBot / headH) / 2;
      if (hTop > 0) Vhead += VheadFull * Math.pow(hTop / headH, 2) * (3 - 2 * hTop / headH) / 2;
    }
  }
  Vtotal = Vcyl + Vhead;
  VtotalFull = VcylFull + VheadFull;
  var pct = VtotalFull > 0 ? Vtotal / VtotalFull * 100 : 0;
  return {Vcyl: Vcyl, Vhead: Vhead, Vtotal: Vtotal, VtotalFull: VtotalFull, pct: pct, Vtotal_L: Vtotal * 1000};
}


var TX = {
  subtitle:{cn:"\u5316\u5b66\u5de5\u7a0b\u8ba1\u7b97\u5e73\u53f0",en:"Chemical Engineering Calculator",bi:"Chemical Engineering Calculator"},
  calcType:{cn:"\u8ba1\u7b97\u76ee\u7684",en:"Calculation Type",bi:"Calc Type"},
  selComp:{cn:"\u9009\u62e9\u7ec4\u5206",en:"Components",bi:"Components \u7ec4\u5206"},
  search:{cn:"\u641c\u7d22\u7ec4\u5206...",en:"Search component...",bi:"Search \u641c\u7d22..."},
  method:{cn:"\u70ed\u529b\u5b66\u65b9\u6cd5",en:"Thermo Method",bi:"Method \u65b9\u6cd5"},
  rec:{cn:"\u63a8\u8350",en:"Rec.",bi:"Rec."},
  cond:{cn:"\u8ba1\u7b97\u6761\u4ef6",en:"Conditions",bi:"Conditions"},
  temp:{cn:"\u6e29\u5ea6",en:"Temp",bi:"Temp T"},
  pres:{cn:"\u538b\u529b",en:"Press",bi:"Press P"},
  calc:{cn:"\u5f00\u59cb\u8ba1\u7b97",en:"Calculate",bi:"Calculate \u8ba1\u7b97"},
  prop:{cn:"\u6027\u8d28",en:"Property",bi:"Property"},
  val:{cn:"\u6570\u503c",en:"Value",bi:"Value"},
  unit:{cn:"\u5355\u4f4d",en:"Unit",bi:"Unit"},
  copy:{cn:"\u590d\u5236CSV",en:"CSV",bi:"CSV"},
  copied:{cn:"\u5df2\u590d\u5236!",en:"Copied!",bi:"OK!"},
  psatCurve:{cn:"\u9971\u548c\u84b8\u6c14\u538b\u66f2\u7ebf",en:"Vapor Pressure Curve",bi:"Psat Curve"},
  empty1:{cn:"\u8f93\u5165\u6e29\u5ea6\u548c\u538b\u529b\uff0c\u8ba1\u7b97\u6c34\u84b8\u6c14\u6027\u8d28",en:"Enter T and P for steam properties",bi:"Enter T, P for steam"},
  empty2:{cn:"\u9009\u62e9\u7ec4\u5206\u5e76\u8f93\u5165\u6761\u4ef6",en:"Select components and conditions",bi:"Select components"},
  iapwsInfo:{cn:"\u57fa\u4e8eIAPWS-IF97\u56fd\u9645\u6807\u51c6",en:"Based on IAPWS-IF97 standard",bi:"IAPWS-IF97 standard"},
  prDesc:{cn:"IGCC/\u5de5\u4e1a\u6c14\u4f53",en:"IGCC & industrial gases",bi:"IGCC gases"},
  srkDesc:{cn:"\u8f7b\u70c3\u4f53\u7cfb",en:"Light hydrocarbons",bi:"Light HCs"},
  ifDesc:{cn:"\u6c34/\u84b8\u6c7d\u6807\u51c6",en:"Water/steam",bi:"Water"},
  idDesc:{cn:"\u4f4e\u538b\u53c2\u8003",en:"Low P ref",bi:"Low P"},
  total:{cn:"\u5408\u8ba1",en:"Total",bi:"Total"},
};
function tx(k, lg) { return TX[k] ? (TX[k][lg] || TX[k].en || k) : k; }
function cl(c, lg) { if (lg === "cn") return c.f+" "+c.cn; if (lg === "en") return c.f+" "+c.en; return c.f+" "+c.cn+"/"+c.en; }

// Fresh color palette: clean light theme with vibrant accents
var C = {
  bg: "#f0f4f8", card: "#ffffff", border: "#e2e8f0", borderH: "#cbd5e1",
  pri: "#0891b2", priL: "#ecfeff", priD: "#0e7490",
  acc: "#f59e0b", accL: "#fffbeb",
  text: "#1e293b", textM: "#475569", textL: "#94a3b8",
  cProp: "#1e40af", cPropL: "#eff6ff",
  cEquip: "#b45309", cEquipL: "#fffbeb",
  cTool: "#047857", cToolL: "#ecfdf5",
  cInst: "#7c3aed", cInstL: "#f5f3ff",
  ok: "#059669", okL: "#ecfdf5",
  err: "#dc2626", errL: "#fef2f2",
  white: "#ffffff",
  shadow: "0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.04)",
  shadowM: "0 4px 12px rgba(0,0,0,0.08)",
  radius: 10,
};

function fmtV(v) {
  if (typeof v !== "number") return String(v);
  if (isNaN(v)) return "\u2014";
  if (Math.abs(v) > 1e6 || (Math.abs(v) < 0.001 && v !== 0)) return v.toExponential(4);
  return v.toPrecision(6);
}

function CompSearch(props) {
  var ref = useRef(null);
  var st = useState(""); var q = st[0]; var setQ = st[1];
  var st2 = useState(false); var open = st2[0]; var setOpen = st2[1];
  var filtered = useMemo(function(){
    var ex = props.exclude || [];
    return COMPS.filter(function(c){ if (ex.indexOf(c.id)>=0) return false; if (!q) return true; var ql=q.toLowerCase(); return c.id.toLowerCase().indexOf(ql)>=0||c.f.toLowerCase().indexOf(ql)>=0||c.cn.indexOf(ql)>=0||c.en.toLowerCase().indexOf(ql)>=0; });
  }, [q, props.exclude]);
  useEffect(function(){ var h=function(e){if(ref.current&&!ref.current.contains(e.target))setOpen(false)}; document.addEventListener("mousedown",h); return function(){document.removeEventListener("mousedown",h)}; },[]);
  return (
    <div ref={ref} style={{position:"relative",flex:1}}>
      <input type="text" value={q} onChange={function(e){setQ(e.target.value);setOpen(true)}} onFocus={function(){setOpen(true)}} placeholder={tx("search",props.lang)}
        style={{width:"100%",padding:"10px 14px",backgroundColor:C.bg,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,outline:"none",boxSizing:"border-box",transition:"border 0.2s"}} />
      {open && filtered.length > 0 ? (
        <div style={{position:"absolute",top:"calc(100% + 4px)",left:0,right:0,backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,maxHeight:240,overflowY:"auto",zIndex:1000,boxShadow:C.shadowM}}>
          {filtered.map(function(c){return (
            <div key={c.id} onClick={function(){props.onSelect(c);setQ("");setOpen(false)}}
              style={{padding:"10px 14px",cursor:"pointer",display:"flex",justifyContent:"space-between",alignItems:"center",borderBottom:"1px solid "+C.border,transition:"background 0.15s"}}
              onMouseEnter={function(e){e.currentTarget.style.backgroundColor=C.priL}}
              onMouseLeave={function(e){e.currentTarget.style.backgroundColor="transparent"}}>
              <span style={{color:C.pri,fontWeight:700,fontFamily:"monospace",fontSize:15}}>{c.f}</span>
              <span style={{color:C.textM,fontSize:12}}>{c.cn} / {c.en}</span>
            </div>
          )})}
        </div>
      ) : null}
    </div>
  );
}

function buildRows(res, lg) {
  if (!res) return [];
  var r = [];
  var addR = function(n,v,u,sep){r.push({n:n,v:v,u:u||"",sep:!!sep})};
  if (res.isVLE) {
    var tp = res.type;
    addR(lg==="en"?"Phase":"\u76f8\u6001", res.ph || "two-phase", "");
    if (tp === "flash") { addR(lg==="en"?"Vapor Fraction":"\u6c14\u76f8\u5206\u7387 VF", res.beta, "--"); addR("T", res.T ? res.T-273.15 : "", "\u00b0C"); addR("P", res.P ? res.P/1e6 : "", "MPa(a)"); }
    if (tp === "bubP") { addR(lg==="en"?"Bubble P":"\u6ce1\u70b9\u538b\u529b", res.P/1e6, "MPa(a)"); addR("T", res.T-273.15, "\u00b0C"); }
    if (tp === "bubT") { addR(lg==="en"?"Bubble T":"\u6ce1\u70b9\u6e29\u5ea6", res.T-273.15, "\u00b0C"); addR("P", res.P/1e6, "MPa(a)"); }
    if (tp === "dewP") { addR(lg==="en"?"Dew P":"\u9732\u70b9\u538b\u529b", res.P/1e6, "MPa(a)"); addR("T", res.T-273.15, "\u00b0C"); }
    if (tp === "dewT") { addR(lg==="en"?"Dew T":"\u9732\u70b9\u6e29\u5ea6", res.T-273.15, "\u00b0C"); addR("P", res.P/1e6, "MPa(a)"); }
    if (tp === "phflash") {
      addR(lg==="en"?"Flash T (result)":"T \u95ea\u84b8\u7ed3\u679c", res.T ? res.T-273.15 : "", "\u00b0C");
      addR(lg==="en"?"Vapor Fraction":"\u6c14\u76f8\u5206\u7387 VF", res.beta, "");
      addR("P", res.P ? res.P/1e6 : "", "MPa(a)");
      addR("H target", res.Htarget, "J/mol");
    }
    if (tp === "psflash") {
      addR(lg==="en"?"Flash T (result)":"T \u95ea\u84b8\u7ed3\u679c", res.T ? res.T-273.15 : "", "\u00b0C");
      addR(lg==="en"?"Vapor Fraction":"\u6c14\u76f8\u5206\u7387 VF", res.beta, "");
      addR("P", res.P ? res.P/1e6 : "", "MPa(a)");
      addR("S target", res.Starget, "J/(mol\u00b7K)");
    }
    if (res.iter) addR(lg==="en"?"Iterations":"\u8fed\u4ee3", res.iter, "");
    if (res.warn) addR(lg==="en"?"Warning":"\u8b66\u544a", res.warn, "");
    if (res.Hmix !== undefined) addR(lg==="en"?"H mix (PR EOS)":"H\u6df7 (PR EOS)", res.Hmix, "J/mol");
    if (res.Smix !== undefined) addR(lg==="en"?"S mix (PR EOS)":"S\u6df7 (PR EOS)", res.Smix, "J/(mol\u00b7K)");
    addR("--- K-values ---","","",true);
    for (var ik=0;ik<(res.clist||[]).length;ik++) addR("K("+res.clist[ik].f+")", res.K?res.K[ik]:"", "");
    if (res.x&&(tp==="flash"||tp==="dewP"||tp==="dewT"||tp==="phflash"||tp==="psflash")) { addR("--- "+(lg==="en"?"Liquid xi":"\u6db2\u76f8 xi")+" ---","","",true); for(var ix=0;ix<res.x.length;ix++) addR(res.clist[ix].f,res.x[ix],""); }
    if (res.y&&(tp==="flash"||tp==="bubP"||tp==="bubT"||tp==="phflash"||tp==="psflash")) { addR("--- "+(lg==="en"?"Vapor yi":"\u6c14\u76f8 yi")+" ---","","",true); for(var iy=0;iy<res.y.length;iy++) addR(res.clist[iy].f,res.y[iy],""); }
    return r;
  }
  if (res.isSteam || res.rg) {
    addR(lg==="en"?"Phase":"\u76f8\u6001",(res.phCN||"")+" / "+(res.ph||""),"");
    if (res.Tresult) addR(lg==="en"?"Temperature (result)":"\u6e29\u5ea6 (\u8ba1\u7b97\u7ed3\u679c)",res.Tresult-273.15,"\u00b0C");
    if (res.vaporFrac!==undefined) addR(lg==="en"?"Vapor Fraction":"\u6c14\u76f8\u5206\u7387 VF",res.vaporFrac,"");
    if (res.Tsat) addR(lg==="en"?"Sat. Temperature":"\u9971\u548c\u6e29\u5ea6 Tsat",res.Tsat-273.15,"\u00b0C");
    if (res.Psat) addR(lg==="en"?"Sat. Pressure":"\u9971\u548c\u538b\u529b Psat",res.Psat/1e6,"MPa(a)");
    if (res.rg===4) {
      addR("--- "+(lg==="en"?"Sat. Liquid":"\u9971\u548c\u6db2\u76f8")+" ---","","",true);
      addR(lg==="en"?"Enthalpy h":"\u7113 h",res.liq.h,"kJ/kg"); addR(lg==="en"?"Entropy s":"\u71b5 s",res.liq.s,"kJ/(kg*K)"); addR(lg==="en"?"Volume v":"\u6bd4\u5bb9 v",res.liq.v,"m3/kg"); addR(lg==="en"?"Density":"\u5bc6\u5ea6",res.liq.rho,"kg/m3"); addR("Cp",res.liq.cp,"kJ/(kg*K)"); addR("Cv",res.liq.cv,"kJ/(kg*K)"); addR(lg==="en"?"Sound Speed":"\u58f0\u901f",res.liq.w,"m/s"); addR(lg==="en"?"Int. Energy u":"\u5185\u80fd u",res.liq.u,"kJ/kg");
      addR("--- "+(lg==="en"?"Sat. Vapor":"\u9971\u548c\u6c14\u76f8")+" ---","","",true);
      addR(lg==="en"?"Enthalpy h":"\u7113 h",res.vap.h,"kJ/kg"); addR(lg==="en"?"Entropy s":"\u71b5 s",res.vap.s,"kJ/(kg*K)"); addR(lg==="en"?"Volume v":"\u6bd4\u5bb9 v",res.vap.v,"m3/kg"); addR(lg==="en"?"Density":"\u5bc6\u5ea6",res.vap.rho,"kg/m3"); addR("Cp",res.vap.cp,"kJ/(kg*K)"); addR("Cv",res.vap.cv,"kJ/(kg*K)"); addR(lg==="en"?"Sound Speed":"\u58f0\u901f",res.vap.w,"m/s"); addR(lg==="en"?"Int. Energy u":"\u5185\u80fd u",res.vap.u,"kJ/kg");
      addR(lg==="en"?"Heat of Vap.":"\u6c7d\u5316\u6f5c\u70ed",res.dhv,"kJ/kg");
    } else {
      addR(lg==="en"?"Enthalpy h":"\u7113 h",res.h,"kJ/kg"); addR(lg==="en"?"Entropy s":"\u71b5 s",res.s,"kJ/(kg*K)"); addR(lg==="en"?"Volume v":"\u6bd4\u5bb9 v",res.v,"m3/kg"); addR(lg==="en"?"Density":"\u5bc6\u5ea6",res.rho,"kg/m3"); addR("Cp",res.cp,"kJ/(kg*K)"); addR("Cv",res.cv,"kJ/(kg*K)"); addR(lg==="en"?"Sound Speed":"\u58f0\u901f",res.w,"m/s"); addR(lg==="en"?"Int. Energy u":"\u5185\u80fd u",res.u,"kJ/kg");
    }
  } else {
    addR(lg==="en"?"Phase":"\u76f8\u6001",(res.phCN||"")+" / "+(res.ph||""),"");
    if (res.vaporFrac!==undefined) addR(lg==="en"?"Vapor Fraction":"\u6c14\u76f8\u5206\u7387",res.vaporFrac,"");
    if (res.Tsat) addR(lg==="en"?"Sat. Temp":"\u9971\u548c\u6e29\u5ea6",res.Tsat-273.15,"\u00b0C");
    if(res.comp&&res.comp.cas)addR("CAS",res.comp.cas,""); if(res.comp&&res.comp.lel)addR("LEL/UEL",res.comp.lel+"/"+res.comp.uel,"vol%"); if(res.comp&&res.comp.fp!==undefined)addR(lg==="en"?"Flash Pt":"\u95ea\u70b9",res.comp.fp,"\u00b0C"); if(res.comp&&res.comp.ait)addR(lg==="en"?"Auto-ignition":"\u81ea\u71c3\u70b9",res.comp.ait,"\u00b0C"); if(res.comp&&res.comp.Tb)addR(lg==="en"?"Boiling Pt":"\u6cb8\u70b9Tb",Math.round((res.comp.Tb-273.15)*100)/100,"\u00b0C"); addR("Z",res.Z,""); addR(lg==="en"?"Molar Vol":"\u6469\u5c14\u4f53\u79ef Vm",res.Vm*1000,"L/mol");
    addR(lg==="en"?"Density":"\u5bc6\u5ea6",res.rho,res.liqDenSrc?"kg/m3 [DIPPR]":"kg/m3");
    addR(lg==="en"?"Enthalpy H":"\u7113 H",res.H/1000,"kJ/mol"); addR(lg==="en"?"Entropy S":"\u71b5 S",res.S,"J/(mol*K)"); addR("Cp",res.Cp,"J/(mol*K)"); addR("Cv",res.gamma>0?res.Cp/res.gamma:0,"J/(mol*K)"); addR("Cp/Cv",res.gamma,""); addR(lg==="en"?"Int. Energy u":"\u5185\u80fd u",res.H-res.Z*8.314*res.T,"J/mol");
    if (res.Psat) addR(lg==="en"?"Sat. Pressure":"\u9971\u548cP",res.Psat/1e6,"MPa(a)");
    if (res.rhoL) addR(lg==="en"?"Sat.Liq.Density":"\u9971\u548c\u6db2\u5bc6\u5ea6",res.rhoL,"kg/m3 [DIPPR]");
    if (res.rhoV) addR(lg==="en"?"Sat.Vap.Density":"\u9971\u548c\u6c14\u5bc6\u5ea6",res.rhoV,"kg/m3");
    addR("--- "+(lg==="en"?"Critical":"Critical \u4e34\u754c")+" ---","","",true);
    addR("Tc",res.comp.Tc,"K"); addR("Pc",res.comp.Pc/1e6,"MPa"); addR("\u03c9",res.comp.w,"");
  }
  return r;
}

function ResTable(props) {
  var cs = useState(false); var copied = cs[0]; var setCopied = cs[1];
  var doCopy = function(){var csv=props.rows.map(function(r){return r.n+","+fmtV(r.v)+","+r.u}).join("\n");navigator.clipboard.writeText("Prop,Value,Unit\n"+csv);setCopied(true);setTimeout(function(){setCopied(false)},2000)};
  return (
    <div style={{backgroundColor:C.white,borderRadius:C.radius,boxShadow:C.shadow,overflow:"hidden",border:"1px solid "+C.border}}>
      {props.title ? <div style={{padding:"12px 16px",display:"flex",justifyContent:"space-between",alignItems:"center",borderBottom:"1px solid "+C.border,background:"linear-gradient(135deg,"+C.priL+","+C.white+")"}}>
        <span style={{color:C.priD,fontWeight:700,fontSize:13}}>{props.title}</span>
        <button onClick={doCopy} style={{padding:"4px 12px",backgroundColor:copied?C.ok:C.bg,border:"1px solid "+(copied?C.ok:C.border),borderRadius:6,color:copied?C.white:C.textM,cursor:"pointer",fontSize:11,fontWeight:600,transition:"all 0.2s"}}>{copied?"OK!":"CSV"}</button>
      </div> : null}
      <table style={{width:"100%",borderCollapse:"collapse"}}>
        <tbody>
          {props.rows.map(function(r,i){return (
            <tr key={i} style={{borderBottom:"1px solid "+(r.sep?"transparent":C.border),backgroundColor:r.sep?C.priL:(i%2===0?C.white:"#f8fafc")}}>
              <td style={{padding:"8px 16px",color:r.sep?C.pri:C.text,fontSize:12,fontWeight:r.sep?700:500}}>{r.n}</td>
              <td style={{padding:"8px 16px",textAlign:"right",color:r.sep?"transparent":C.pri,fontFamily:"monospace",fontWeight:700,fontSize:13}}>{fmtV(r.v)}</td>
              <td style={{padding:"8px 16px",color:C.textL,fontSize:11}}>{r.u}</td>
            </tr>
          )})}
        </tbody>
      </table>
    </div>
  );
}

// Pill button helper
function Pill(props) {
  var active = props.active;
  return (
    <button onClick={props.onClick} disabled={props.disabled}
      style={{padding:"8px 14px",borderRadius:20,border:active?"2px solid "+C.pri:"1.5px solid "+C.border,
        backgroundColor:active?C.pri:C.white,color:active?C.white:C.textM,
        cursor:props.disabled?"not-allowed":"pointer",fontSize:12,fontWeight:600,transition:"all 0.2s",
        boxShadow:active?"0 2px 8px rgba(8,145,178,0.3)":"none",flex:"1 1 auto",textAlign:"center",
        opacity:props.disabled?0.4:1}}>
      {props.children}
    </button>
  );
}

export default function App() {
  // Load component database from API on first render
  var sDataReady = useState(false), dataReady = sDataReady[0], setDataReady = sDataReady[1];
  var sDataErr = useState(false), dataErr = sDataErr[0], setDataErr = sDataErr[1];
  var sLang = useState("cn"), lang = sLang[0], setLang = sLang[1];
  var sMode = useState("pure"), mode = sMode[0], setMode = sMode[1];
  var sGroup = useState("props"), group = sGroup[0], setGroup = sGroup[1];
  var sComps = useState([]), comps = sComps[0], setComps = sComps[1];
  var sMeth = useState("PR"), meth = sMeth[0], setMeth = sMeth[1];
  var sTv = useState("25"), tv = sTv[0], setTv = sTv[1];
  var sTu = useState("C"), tu = sTu[0], setTu0 = sTu[1];
  var sPv = useState("0.1013"), pv = sPv[0], setPv = sPv[1];
  var sPu = useState("MPa"), pu = sPu[0], setPu0 = sPu[1];

  // Smart unit switch: auto-convert displayed value when unit changes
  var setTu = function(newU) {
    var oldU = TU.find(function(u){return u.id===tu});
    var newUO = TU.find(function(u){return u.id===newU});
    var val = parseFloat(tv);
    if (oldU && newUO && !isNaN(val)) {
      var K = oldU.toK(val);
      var converted = newUO.fromK(K);
      setTv(String(Math.round(converted*1e6)/1e6));
    }
    setTu0(newU);
  };
  var setPu = function(newU) {
    var oldU = PU.find(function(u){return u.id===pu});
    var newUO = PU.find(function(u){return u.id===newU});
    var val = parseFloat(pv);
    if (oldU && newUO && !isNaN(val)) {
      var Pa = oldU.toP(val);
      var converted = newUO.fromP(Pa);
      setPv(String(Math.round(converted*1e8)/1e8));
    }
    setPu0(newU);
  };
  var sVf = useState(""), vf = sVf[0], setVf = sVf[1];
  var sInMode = useState("pt"), inMode = sInMode[0], setInMode = sInMode[1];
  var sHval = useState("3400"), hval = sHval[0], setHval = sHval[1]; // kJ/kg for PH flash
  var sSval = useState("7.0"), sval = sSval[0], setSval = sSval[1]; // kJ/(kg*K) for PS flash
  // COP inputs
  var sCop = useState({tevap:"-10",tcond:"45",tsh:"5",tsc:"3",eta:"0.75"}), copIn = sCop[0];
  var setCop = function(k,v){sEqRes[1](null);sCop[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  // LMTD inputs
  var sLmtd = useState({thi:"120",tho:"80",tci:"20",tco:"60",flow:"counter"}), lmtdIn = sLmtd[0];
  var setLmtd = function(k,v){sEqRes[1](null);sLmtd[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sVleType = useState("flash"), vleType = sVleType[0], setVleType = sVleType[1];
  var sBasis = useState("mol"), basis = sBasis[0], setBasis = sBasis[1];
  var sDiag = useState(null), diagData = sDiag[0], setDiagData = sDiag[1];
  // Equipment inputs
  var sEqRes = useState(null), eqRes = sEqRes[0], setEqRes = sEqRes[1];
  var sEqIn = useState({p1:"0.1013",p2:"0.5",t1:"25",eta:"0.75",gamma:"1.4",mw:"28.97",flow:"100",
    dia:"0.1",len:"100",nst:"1",tcool:"35",rough:"0.000046",vel:"3",rho:"1000",mu:"0.001",qm:"",qv:""}), eqIn = sEqIn[0];
  var setEq = function(k,v){setEqRes(null);sEqIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  // Unit converter
  var sUcCat = useState("temp"), ucCat = sUcCat[0], setUcCat = sUcCat[1];
  var sUcFrom = useState("C"), ucFrom = sUcFrom[0], setUcFrom = sUcFrom[1];
  var sUcTo = useState("K"), ucTo = sUcTo[0], setUcTo = sUcTo[1];
  var sUcVal = useState("100"), ucVal = sUcVal[0], setUcVal = sUcVal[1];
  var sInfo = useState(false), showInfo = sInfo[0], setShowInfo = sInfo[1];
  var sPwIdx = useState(0), pwIdx = sPwIdx[0], setPwIdx = sPwIdx[1];
  var sPwMat = useState("A106B"), pwMat = sPwMat[0], setPwMat = sPwMat[1];
  var sPwEw = useState("1.0"), pwEw = sPwEw[0], setPwEw = sPwEw[1];
  var sPwYc = useState("0.4"), pwYc = sPwYc[0], setPwYc = sPwYc[1];
  var sPwCc = useState("1.5"), pwCc = sPwCc[0], setPwCc = sPwCc[1];
  var sSvType = useState("gas"), svType = sSvType[0], setSvType = sSvType[1];
  var sSvIn = useState({w:"5000",t1:"200",z:"1",mw:"28.97",gamma:"1.4",p1:"1000",pb:"101.3",q:"50",sg:"1.0",kd:"0.975"}), svIn = sSvIn[0];
  var setSvI = function(k,v){setEqRes(null);sSvIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sOpIn = useState({D:"100",d:"50",dP:"5000",rho:"1.2",mu:"0.000018",tap:"flange"}), opIn = sOpIn[0];
  var setOp = function(k,v){setEqRes(null);sOpIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sGfIn = useState({am:"1000",t:"25",p:"0.5",nm:"",sm:"",scfm:"",mmscfd:"",activeField:"am"}), gfIn = sGfIn[0];
  var sCvType = useState("liq"), cvType = sCvType[0], setCvType = sCvType[1];
  var sCvIn = useState({q:"50",g:"1.0",p1:"500",p2:"200",w:"5000",t1:"200",mw:"28.97",gamma:"1.4",z:"1"}), cvIn = sCvIn[0];
  var setCvI = function(k,v){setEqRes(null);sCvIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sVesIn = useState({d:"2",len:"5",h:"1",head:"ellip",orient:"horiz"}), vesIn = sVesIn[0];
  var sInsIn = useState({d:"100",tp:"200",ta:"25",thick:"50",k:"0.04",wind:"3",emis:"0.9",mat:"rockwool"}), insIn = sInsIn[0];
  var setInsI = function(k,v){setEqRes(null);sInsIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sPhaseSel = useState("H2O"), phaseSel = sPhaseSel[0], setPhaseSel = sPhaseSel[1];
  var sMaIn = useState({pv:"50",lo:"0",hi:"100",ma:"12",dir:"pv2ma"}), maIn = sMaIn[0];
  var setMaI = function(k,v){sMaIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sTcType = useState("K"), tcType = sTcType[0], setTcType = sTcType[1];
  var sTcVal = useState("100"), tcVal = sTcVal[0], setTcVal = sTcVal[1];
  var sTcDir = useState("t2mv"), tcDir = sTcDir[0], setTcDir = sTcDir[1];
  var sHxIn = useState({q:"500",u:"1500",lmtd:"30"}), hxIn = sHxIn[0];
  var setHxI = function(k,v){sHxIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sNpshIn = useState({ps:"101.325",pv:"2.339",hs:"3",hf:"0.5",rho:"1000"}), npshIn = sNpshIn[0];
  var setNpshI = function(k,v){sNpshIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sHumIn = useState({t:"25",rh:"50",patm:"101.325"}), humIn = sHumIn[0];
  var setHumI = function(k,v){sHumIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sBd = useState({type:"gas",std:"GB",w:"100",t1:"50",z:"1.0",mw:"28.97",gamma:"1.4",p0:"500",kd:"0.62",q:"100",sg:"1.0",p1:"500",p2:"101.3",kw:"1.0",kv:"1.0"}), bdIn = sBd[0];
  var setBd = function(k,v){sBd[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var sMix = useState({t:"100",p:"0.5",phase:"both"}), mixIn = sMix[0];
  var setMixI = function(k,v){sMix[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var setVesI = function(k,v){setEqRes(null);sVesIn[1](function(p){var o={};for(var x in p)o[x]=p[x];o[k]=v;return o})};
  var setGf = function(k,v,bulk){sGfIn[1](function(p){var o={};for(var x in p)o[x]=p[x];if(bulk){for(var bk in bulk)o[bk]=bulk[bk];}else{o[k]=v;}return o})};
  var sRes = useState(null), res = sRes[0], setRes = sRes[1];
  var sErr = useState(null), err = sErr[0], setErr = sErr[1];
  var sSat = useState(null), satC = sSat[0], setSatC = sSat[1];

  var addComp = useCallback(function(c){setComps(function(p){if(p.find(function(x){return x.id===c.id}))return p;var nc={};for(var k in c)nc[k]=c[k];nc.mf="";return p.concat([nc])});},[]);

  // Persistence: save inputs to storage (Artifact API or localStorage fallback)
  var _stgSet = function(k,v){try{if(typeof window.storage!=="undefined"&&window.storage.set){window.storage.set(k,v);}else{localStorage.setItem(k,v);}}catch(e){}};
  var _stgGet = function(k,cb){try{if(typeof window.storage!=="undefined"&&window.storage.get){window.storage.get(k).then(function(r){cb(r&&r.value?r.value:null)}).catch(function(){});}else{cb(localStorage.getItem(k));}}catch(e){}};
  useEffect(function(){
    _stgSet("ctp_inputs", JSON.stringify({tv:tv,tu:tu,pv:pv,pu:pu,lang:lang,mode:mode,group:group,meth:meth}));
  },[tv,tu,pv,pu,lang,mode,group,meth]);
  // Load component data from Worker API on mount
  useEffect(function(){
    loadChemData().then(function(){
      setDataReady(true);
    }).catch(function(){
      setDataErr(true);
      setDataReady(true); // still show UI, just empty components
    });
  },[]);

  // Load on mount
  useEffect(function(){
    _stgGet("ctp_inputs", function(val){
      if(!val) return;
      try{var d=JSON.parse(val);if(d.tv)setTv(d.tv);if(d.tu)setTu0(d.tu);if(d.pv)setPv(d.pv);if(d.pu)setPu0(d.pu);if(d.lang)setLang(d.lang);if(d.mode)setMode(d.mode);if(d.group)setGroup(d.group);if(d.meth)setMeth(d.meth);}catch(e){}
    });
  },[]);
  var removeComp = useCallback(function(id){setComps(function(p){return p.filter(function(x){return x.id!==id})});},[]);
  var updMf = useCallback(function(id,v){setComps(function(p){return p.map(function(x){return x.id===id?Object.assign({},x,{mf:v}):x})});},[]);
  var totalMf = useMemo(function(){return comps.reduce(function(s,c){return s+(parseFloat(c.mf)||0)},0)},[comps]);

  var doCalc = useCallback(function(){
    setErr(null);setRes(null);setSatC(null);setDiagData(null);
    var tObj=TU.find(function(u){return u.id===tu}), pObj=PU.find(function(u){return u.id===pu});
    // PH flash (steam only)
    if (inMode==="ph" && (mode==="steam"||(comps.length===1&&comps[0]&&comps[0].id==="H2O"))) {
      var PPaPH=pObj.toP(parseFloat(pv)); var Hv=parseFloat(hval);
      if(isNaN(PPaPH)||PPaPH<=0){setErr("Invalid P");return;}
      if(isNaN(Hv)){setErr("Invalid H (kJ/kg)");return;}
      var phR=steamPHflash(PPaPH,Hv);
      if(phR&&phR.props){var sr=phR.props;sr.isSteam=true;sr.inputMode="PH";sr.Tresult=phR.T;setRes(sr);}
      else{setErr("PH flash failed");}
      return;
    }
    // PS flash (steam only)
    if (inMode==="ps" && (mode==="steam"||(comps.length===1&&comps[0]&&comps[0].id==="H2O"))) {
      var PPaPS=pObj.toP(parseFloat(pv)); var Sv=parseFloat(sval);
      if(isNaN(PPaPS)||PPaPS<=0){setErr("Invalid P");return;}
      if(isNaN(Sv)){setErr("Invalid S (kJ/(kg*K))");return;}
      var psR=steamPSflash(PPaPS,Sv);
      if(psR&&psR.props){var sr2=psR.props;sr2.isSteam=true;sr2.inputMode="PS";sr2.Tresult=psR.T;setRes(sr2);}
      else{setErr("PS flash failed");}
      return;
    }
    if ((inMode==="pvf"||inMode==="tvf")&&(mode==="pure"||mode==="steam")) {
      var vfVal=parseFloat(vf); if(isNaN(vfVal)||vfVal<0||vfVal>1){setErr("VF: 0~1");return;}
      if (inMode==="pvf") {
        var PPa2=pObj.toP(parseFloat(pv)); if(isNaN(PPa2)||PPa2<=0){setErr("Invalid P");return;}
        if (mode==="steam"||(comps.length===1&&comps[0]&&comps[0].id==="H2O")) {
          var Tsat=iapws97_Tsat(PPa2); if(!Tsat){setErr("No Tsat");return;}
          var sr=calcSteam(Tsat,PPa2); if(sr.err){setErr(sr.err);return;}
          sr.vaporFrac=vfVal;sr.Tsat=Tsat; setRes(Object.assign({isSteam:true},sr)); return;
        }
        var c0=comps[0]; if(!c0){setErr("Select component");return;}
        var Tlo=0.35*c0.Tc,Thi=c0.Tc*0.999;
        for(var bi=0;bi<100;bi++){var Tm=(Tlo+Thi)/2;var Pm=calcPsat(Tm,c0,meth);if(!Pm){Tlo=Tm;continue;}if(Pm<PPa2)Tlo=Tm;else Thi=Tm;if(Math.abs(Pm-PPa2)/PPa2<1e-8)break;}
        var Ts2=(Tlo+Thi)/2;var r2=calcPure(Ts2,PPa2,c0,meth);r2.ph="two-phase";r2.phCN="\u4e24\u76f8\u533a";r2.vaporFrac=vfVal;r2.Tsat=Ts2;setRes(r2);return;
      }
      if (inMode==="tvf") {
        var TK2=tObj.toK(parseFloat(tv)); if(isNaN(TK2)||TK2<=0){setErr("Invalid T");return;}
        if (mode==="steam"||(comps.length===1&&comps[0]&&comps[0].id==="H2O")) {
          var Ps3=steamPsat(TK2); if(!Ps3){setErr("T out of range");return;}
          var sr3=calcSteam(TK2,Ps3); if(sr3.err){setErr(sr3.err);return;}
          sr3.vaporFrac=vfVal;sr3.Psat=Ps3;sr3.Tsat=TK2;setRes(Object.assign({isSteam:true},sr3));return;
        }
        var c1=comps[0]; if(!c1){setErr("Select component");return;}
        var Ps4=calcPsat(TK2,c1,meth); if(!Ps4){setErr("No Psat");return;}
        var r3=calcPure(TK2,Ps4,c1,meth);r3.ph="two-phase";r3.phCN="\u4e24\u76f8\u533a";r3.vaporFrac=vfVal;r3.Psat=Ps4;r3.Tsat=TK2;setRes(r3);return;
      }
    }
    var TK=tObj.toK(parseFloat(tv)), PPa=pObj.toP(parseFloat(pv));
    // Skip validation for inputs not needed by current mode
    var needT=true, needP=true;
    if(mode==="diag"){needT=vleType==="pxy";needP=vleType==="txy"||vleType==="flash";}
    if(mode==="vle"){needT=vleType==="flash"||vleType==="bubP"||vleType==="dewP";needP=vleType==="flash"||vleType==="bubT"||vleType==="dewT";}
    if(needT&&(isNaN(TK)||TK<=0)){setErr("Invalid T");return;}
    if(needP&&(isNaN(PPa)||PPa<=0)){setErr("Invalid P");return;}
    if(mode==="steam"){var sr4=calcSteam(TK,PPa);if(sr4.err){setErr(sr4.err);return;}setRes(Object.assign({isSteam:true},sr4));return;}
    if(comps.length===0){setErr("Select component");return;}
    if (mode==="vle") {
      if(comps.length<2){setErr("VLE: 2+ components");return;}
      var z=comps.map(function(c){return parseFloat(c.mf)||0});var zS=z.reduce(function(a,b){return a+b},0);
      var autoNorm = Math.abs(zS-1)>0.005 && zS>0;
      if(zS<=0){setErr("All compositions are zero");return;}
      for(var iz=0;iz<z.length;iz++)z[iz]/=zS;
      if(basis==="wt"){var nf=[],ns=0;for(var iw=0;iw<z.length;iw++){var ni=z[iw]/comps[iw].MW;nf.push(ni);ns+=ni}for(var jw=0;jw<z.length;jw++)z[jw]=nf[jw]/ns}
      var clist=comps.map(function(c){return{id:c.id,Tc:c.Tc,Pc:c.Pc,w:c.w,MW:c.MW,f:c.f,cn:c.cn,en:c.en}});
      var vr;
      if(vleType==="flash"){
        vr=tpFlash(TK,PPa,z,clist);vr.type="flash";vr.T=TK;vr.P=PPa;
        // Q3: compute mixture H and S for TP flash
        try{var HL=mixEnthalpy(TK,PPa,vr.x||z,clist,"liquid"),HV=mixEnthalpy(TK,PPa,vr.y||z,clist,"vapor");vr.Hmix=(1-vr.beta)*HL+vr.beta*HV;}catch(e){}
        try{var SL=mixEntropy(TK,PPa,vr.x||z,clist,"liquid"),SV=mixEntropy(TK,PPa,vr.y||z,clist,"vapor");vr.Smix=(1-vr.beta)*SL+vr.beta*SV;}catch(e){}
      }
      else if(vleType==="bubP"){vr=bubbleP(TK,z,clist);vr.type="bubP";vr.T=TK;}
      else if(vleType==="bubT"){vr=bubbleT(PPa,z,clist);vr.type="bubT";vr.P=PPa;}
      else if(vleType==="dewP"){vr=dewP(TK,z,clist);vr.type="dewP";vr.T=TK;}
      else if(vleType==="dewT"){vr=dewT(PPa,z,clist);vr.type="dewT";vr.P=PPa;}
      else if(vleType==="phflash"){
        var Htgt=parseFloat(hval); if(isNaN(Htgt)){setErr("Enter H (J/mol)");return;}
        vr=vlePhFlash(PPa,Htgt,z,clist);vr.type="phflash";vr.P=PPa;vr.Htarget=Htgt;
      }
      else if(vleType==="psflash"){
        var Stgt=parseFloat(sval); if(isNaN(Stgt)){setErr("Enter S (J/(mol\u00b7K))");return;}
        vr=vlePsFlash(PPa,Stgt,z,clist);vr.type="psflash";vr.P=PPa;vr.Starget=Stgt;
      }
      vr.isVLE=true;vr.clist=clist;vr.z=z;vr.autoNorm=autoNorm;setRes(vr);return;
    }
    if (mode==="diag") {
      if(comps.length!==2){setErr("Diagram: 2 components");return;}
      var d1={id:comps[0].id,Tc:comps[0].Tc,Pc:comps[0].Pc,w:comps[0].w,MW:comps[0].MW,f:comps[0].f};
      var d2={id:comps[1].id,Tc:comps[1].Tc,Pc:comps[1].Pc,w:comps[1].w,MW:comps[1].MW,f:comps[1].f};
      if(vleType==="txy"||vleType==="flash"){var dd=genTxy(PPa,d1,d2,30);setDiagData({type:"txy",bub:dd.bub,dew:dd.dew,P:PPa,c1:comps[0],c2:comps[1],supercrit:dd.supercrit,Ts1:dd.Ts1,Ts2:dd.Ts2});}
      else{var dd2=genPxy(TK,d1,d2,30);setDiagData({type:"pxy",bub:dd2.bub,dew:dd2.dew,T:TK,c1:comps[0],c2:comps[1],supercrit:dd2.supercrit,Ps1:dd2.Ps1,Ps2:dd2.Ps2});}
      return;
    }
    var c0b=comps[0];
    if(c0b.id==="H2O"&&comps.length===1&&meth!=="SRK"){var sr5=calcSteam(TK,PPa);if(!sr5.err){setRes(Object.assign({isSteam:true,comp:c0b},sr5));return;}}
    var rb=calcPure(TK,PPa,c0b,meth);if(rb.err){setErr(rb.err);return;}setRes(rb);
    if(TK<c0b.Tc){var dd3=[];var curTC=TK-273.15;var tMn=Math.max(curTC-80,c0b.Tc*0.35-273.15);var tMx=Math.min(curTC+80,c0b.Tc*0.999-273.15);if(tMx-tMn<20){tMn=c0b.Tc*0.35-273.15;tMx=c0b.Tc*0.999-273.15;}for(var ii=0;ii<=60;ii++){var tC=tMn+(tMx-tMn)*ii/60;var tK=tC+273.15;if(tK<c0b.Tc*0.35||tK>c0b.Tc*0.999)continue;var ps=calcPsat(tK,c0b,meth);if(ps&&ps>0)dd3.push({T:Math.round(tC*10)/10,Psat:Math.round(ps/1e6*1e6)/1e6});}setSatC(dd3);}
  },[mode,comps,meth,tv,tu,pv,pu,vf,inMode,vleType,hval,sval,basis]);

  var rows = useMemo(function(){return buildRows(res,lang)},[res,lang]);
  var modes=[{id:"pure",g:"props",cn:"\u7eaf\u7ec4\u5206",en:"Pure"},{id:"steam",g:"props",cn:"\u84b8\u6c14",en:"Steam"},{id:"vle",g:"props",cn:"VLE",en:"VLE"},{id:"diag",g:"props",cn:"\u76f8\u56fe",en:"Diagram"},{id:"comp",g:"equip",cn:"\u538b\u7f29\u673a",en:"Compressor"},{id:"expd",g:"equip",cn:"\u81a8\u80c0\u673a",en:"Expander"},{id:"pump",g:"equip",cn:"\u6cf5",en:"Pump"},{id:"pipe",g:"equip",cn:"\u7ba1\u8def",en:"Pipe"},{id:"cop",g:"equip",cn:"\u5236\u51b7\u5faa\u73af",en:"Refrig. COP"},{id:"lmtd",g:"equip",cn:"LMTD",en:"LMTD"},{id:"unit",g:"tools",cn:"\u5355\u4f4d\u6362\u7b97",en:"Unit Conv."},{id:"pwall",g:"tools",cn:"\u7ba1\u9053\u58c1\u539a",en:"Pipe Wall"},{id:"svsize",g:"inst",cn:"\u5b89\u5168\u9600",en:"Safety Valve"},{id:"oplate",g:"inst",cn:"\u5b54\u677f",en:"Orifice"},{id:"gasflow",g:"tools",cn:"\u6807\u51c6\u6d41\u91cf",en:"Std Gas Flow"},{id:"mwcalc",g:"tools",cn:"\u5206\u5b50\u91cf",en:"MW Calc"},{id:"cvsize",g:"inst",cn:"\u9600\u95e8Cv",en:"Valve Cv"},{id:"heatval",g:"tools",cn:"\u70ed\u503c",en:"Heat Value"},{id:"vessel",g:"tools",cn:"\u5bb9\u5668\u6db2\u4f4d",en:"Vessel Vol."},{id:"hxarea",g:"equip",cn:"\u6362\u70ed\u5668",en:"HX Area"},{id:"npsh",g:"equip",cn:"NPSH",en:"NPSH"},{id:"ma420",g:"inst",cn:"4-20mA",en:"4-20mA"},{id:"tcconv",g:"inst",cn:"\u70ed\u7535\u5076",en:"T/C"},{id:"humid",g:"tools",cn:"\u6e7f\u5ea6\u9732\u70b9",en:"Humidity"},{id:"b3610",g:"tools",cn:"\u7ba1\u89c4B36",en:"Pipe Sch."},{id:"pipevel",g:"tools",cn:"\u6d41\u901f\u8868",en:"Velocity"},{id:"insul",g:"equip",cn:"\u4fdd\u6e29\u6563\u70ed",en:"Insulation"},{id:"phasediag",g:"props",cn:"\u76f8\u56fe(PT)",en:"Phase(PT)"},{id:"mixprop",g:"props",cn:"\u6df7\u5408\u7269\u6027",en:"Mix Props"},{id:"burstdisc",g:"inst",cn:"\u7206\u7834\u7247",en:"Burst Disc"}];
  var groups=[{id:"props",cn:"\u7269\u6027\u8ba1\u7b97",en:"Properties",bi:"Properties \u7269\u6027"},{id:"equip",cn:"\u8bbe\u5907\u8ba1\u7b97",en:"Equipment",bi:"Equipment \u8bbe\u5907"},{id:"inst",cn:"\u4eea\u8868\u81ea\u63a7",en:"Instrument"},{id:"tools",cn:"\u5de5\u5177",en:"Tools",bi:"Tools \u5de5\u5177"}];
  var meths=[{id:"PR",nm:"Peng-Robinson",ds:tx("prDesc",lang)},{id:"SRK",nm:"SRK",ds:tx("srkDesc",lang)},{id:"IAPWS-IF97",nm:"IAPWS-IF97",ds:tx("ifDesc",lang)},{id:"Ideal",nm:"Ideal Gas",ds:tx("idDesc",lang)}];
  var recM=mode==="steam"?"IAPWS-IF97":(comps.length===1&&comps[0]&&comps[0].id==="H2O")?"IAPWS-IF97":"PR";

  var inputStyle = {flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
  var selStyle = {padding:"10px 8px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:12,cursor:"pointer",outline:"none"};
  var cardStyle = {backgroundColor:C.white,borderRadius:C.radius+2,padding:16,marginBottom:12,boxShadow:C.shadow,border:"1px solid "+C.border};

  // Determine which inputs to show
  var showT=true,showP=true,showVF=false,showH=false,showS=false;
  if(mode==="pure"||mode==="steam"){showT=inMode!=="pvf"&&inMode!=="ph"&&inMode!=="ps";showP=inMode!=="tvf";showVF=inMode==="pvf"||inMode==="tvf";showH=inMode==="ph";showS=inMode==="ps";}
  else if(mode==="vle"){
    var isVlePH=vleType==="phflash", isVlePS=vleType==="psflash";
    showT=!isVlePH&&!isVlePS&&(vleType==="flash"||vleType==="bubP"||vleType==="dewP");
    showP=vleType==="flash"||vleType==="bubT"||vleType==="dewT"||isVlePH||isVlePS;
    showH=isVlePH; showS=isVlePS;
  }
  else if(mode==="diag"){showT=vleType==="pxy";showP=vleType==="txy"||vleType==="flash";}

  // Show loading screen while fetching component data from Worker API
  if (!dataReady) return (
    <div style={{minHeight:"100vh",backgroundColor:"#0f172a",display:"flex",alignItems:"center",justifyContent:"center",flexDirection:"column",gap:16}}>
      <div style={{fontSize:28,fontWeight:700,color:"#0891b2",letterSpacing:1}}>ChemCalc</div>
      <div style={{color:"#94a3b8",fontSize:13}}>{lang==="en"?"Loading component database...":"正在加载组分数据库..."}</div>
      <div style={{width:40,height:40,border:"3px solid #1e293b",borderTop:"3px solid #0891b2",borderRadius:"50%",animation:"spin 0.8s linear infinite"}} />
      <style>{"@keyframes spin{to{transform:rotate(360deg)}}"}</style>
      {dataErr?<div style={{color:"#ef4444",fontSize:11,marginTop:8}}>{lang==="en"?"Network error - check internet connection":"网络错误，请检查网络连接"}</div>:null}
    </div>
  );

  return (
    <div style={{minHeight:"100vh",backgroundColor:C.bg,color:C.text,fontFamily:"-apple-system,sans-serif"}}>
      {/* Header */}
      <div style={{background:"linear-gradient(135deg, #0f172a 0%, #1e293b 50%, #334155 100%)",padding:"18px 20px",display:"flex",justifyContent:"space-between",alignItems:"center",flexWrap:"wrap",gap:8,boxShadow:"0 2px 12px rgba(8,145,178,0.25)"}}>
        <div>
          <div style={{fontSize:22,fontWeight:700,color:"#ffffff",letterSpacing:0.5}}>ChemCalc</div>
          <div style={{fontSize:11,color:"rgba(255,255,255,0.7)",marginTop:1}}>{tx("subtitle",lang)} | PR / SRK / IAPWS-IF97</div>
        </div>
        <div style={{display:"flex",gap:4,alignItems:"center"}}>
          {[{id:"cn",lb:"\u4e2d\u6587"},{id:"en",lb:"EN"}].map(function(b){return (
            <button key={b.id} onClick={function(){setLang(b.id)}} style={{padding:"5px 12px",backgroundColor:lang===b.id?"rgba(255,255,255,0.25)":"rgba(255,255,255,0.1)",border:"1px solid rgba(255,255,255,"+(lang===b.id?"0.5":"0.2")+")",borderRadius:16,color:"#fff",cursor:"pointer",fontSize:11,fontWeight:lang===b.id?700:400,transition:"all 0.2s"}}>{b.lb}</button>
          )})}
          <button onClick={function(){setShowInfo(true)}} style={{padding:"5px 12px",backgroundColor:"rgba(255,255,255,0.1)",border:"1px solid rgba(255,255,255,0.2)",borderRadius:16,color:"#fff",cursor:"pointer",fontSize:11,fontWeight:600}}>{lang==="en"?"About":"\u5173\u4e8e"}</button>
        </div>
      </div>

      <div style={{maxWidth:1280,margin:"0 auto",padding:"16px 12px"}}>
        {/* Level 1: Group tabs */}
        <div style={{display:"flex",gap:4,marginBottom:8}}>
          {groups.map(function(g){return (
            <button key={g.id} onClick={function(){setGroup(g.id);var first=modes.filter(function(m){return m.g===g.id})[0];if(first)setMode(first.id);setRes(null);setErr(null);setDiagData(null);setEqRes(null)}}
              style={{flex:1,padding:"10px 8px",borderRadius:8,border:group===g.id?"2px solid "+(g.id==="props"?C.cProp:g.id==="equip"?C.cEquip:g.id==="inst"?C.cInst:C.cTool):"1.5px solid "+C.border,backgroundColor:group===g.id?(g.id==="props"?C.cProp:g.id==="equip"?C.cEquip:g.id==="inst"?C.cInst:C.cTool):C.white,color:group===g.id?"#fff":C.textM,cursor:"pointer",fontSize:13,fontWeight:700,transition:"all 0.2s"}}>
              {g[lang]||g.en}
            </button>
          )})}
        </div>
        {/* Level 2: Sub-mode pills */}
        {/* Category accent color */}
        <div style={{display:"flex",gap:4,marginBottom:14,flexWrap:"wrap"}}>
          {modes.filter(function(m){return m.g===group}).map(function(m){return (
            <Pill key={m.id} active={mode===m.id} onClick={function(){setMode(m.id);if(m.id==="steam"){setMeth("IAPWS-IF97")}else{if(inMode==="ph"||inMode==="ps")setInMode("pt");}setRes(null);setErr(null);setDiagData(null);setEqRes(null);if(m.g!=="props")setComps([])}}>
              {m[lang]||m.en}
            </Pill>
          )})}
        </div>

        <div style={{display:"flex",gap:14,flexWrap:"wrap",alignItems:"flex-start"}}>
          {/* LEFT PANEL - only for Properties group */}
          {group==="props"?(<div style={{flex:"1 1 380px",minWidth:300,maxWidth:480}}>
            {/* Components */}
            {mode!=="steam"?(
              <div style={cardStyle}>
                <div style={{fontSize:11,color:C.textL,marginBottom:8,fontWeight:600,textTransform:"uppercase",letterSpacing:1}}>{tx("selComp",lang)}</div>
                {comps.map(function(c,idx){return (
                  <div key={c.id} style={{display:"flex",alignItems:"center",gap:6,marginBottom:6,padding:"8px 10px",backgroundColor:C.bg,borderRadius:8,border:"1px solid "+C.border}}>
                    <span style={{color:C.pri,fontWeight:700,fontFamily:"monospace",minWidth:55,fontSize:14}}>{c.f}</span>
                    <span style={{color:C.textM,fontSize:11,flex:1}}>{cl(c,lang)}</span>
                    <span style={{color:C.textL,fontSize:9,minWidth:50,textAlign:"right"}}>{"Tc="+(c.Tc>200?(c.Tc-273.15).toFixed(0)+"\u00b0C":(c.Tc-273.15).toFixed(0)+"\u00b0C")}</span>
                    {comps.length>1?<input type="number" step="0.01" min="0" max="1" placeholder="xi" value={c.mf} onChange={function(e){updMf(c.id,e.target.value)}} style={{width:58,padding:"4px 6px",backgroundColor:C.white,border:"1px solid "+C.border,borderRadius:6,color:C.text,fontSize:12,textAlign:"center",fontFamily:"monospace"}} />:null}
                    <button onClick={function(){removeComp(c.id)}} style={{background:"none",border:"none",color:C.err,cursor:"pointer",fontSize:16,padding:"0 4px",fontWeight:700}}>x</button>
                  </div>
                )})}
                {(mode==="pure"&&comps.length>=1)||(mode==="diag"&&comps.length>=2)?null:<CompSearch onSelect={addComp} lang={lang} exclude={comps.map(function(c){return c.id})} />}
                {comps.length>1?(<div>
                  <div style={{display:"flex",gap:4,marginTop:6,marginBottom:4}}>
                    <Pill active={basis==="mol"} onClick={function(){setBasis("mol")}}>mol%</Pill>
                    <Pill active={basis==="wt"} onClick={function(){setBasis("wt")}}>wt%</Pill>
                  </div>
                  <div style={{padding:"5px 10px",borderRadius:6,backgroundColor:Math.abs(totalMf-1)<0.005?C.okL:C.accL,color:Math.abs(totalMf-1)<0.005?C.ok:C.acc,fontSize:11,display:"flex",justifyContent:"space-between",fontWeight:600}}>
                  <span>{tx("total",lang)} ({basis})</span><span style={{fontFamily:"monospace"}}>{"= "+totalMf.toFixed(4)}</span>
                  </div>
                </div>):null}
              </div>
            ):null}

            {/* Method */}
            {mode!=="steam"?(
              <div style={cardStyle}>
                <div style={{fontSize:11,color:C.textL,marginBottom:8,fontWeight:600,textTransform:"uppercase",letterSpacing:1}}>{tx("method",lang)}</div>
                <div style={{display:"flex",gap:6,flexWrap:"wrap"}}>
                  {meths.map(function(m){return (
                    <div key={m.id} onClick={function(){setMeth(m.id)}} style={{flex:"1 1 auto",padding:"8px 12px",borderRadius:8,cursor:"pointer",backgroundColor:meth===m.id?C.priL:C.bg,border:meth===m.id?"1.5px solid "+C.pri:"1.5px solid transparent",transition:"all 0.2s",textAlign:"center"}}>
                      <div style={{fontSize:12,fontWeight:700,color:meth===m.id?C.pri:C.text}}>{m.nm}</div>
                      <div style={{fontSize:9,color:C.textL,marginTop:2}}>{m.ds} {recM===m.id?" *":""}</div>
                    </div>
                  )})}
                </div>
              </div>
            ):null}

            {/* VLE/Diagram type selector */}
            {(mode==="vle"||mode==="diag")?(
              <div style={cardStyle}>
                <div style={{fontSize:11,color:C.textL,marginBottom:8,fontWeight:600,textTransform:"uppercase",letterSpacing:1}}>
                  {mode==="vle"?(lang==="en"?"VLE Type":"相平衡类型"):(lang==="en"?"Diagram Type":"相图类型")}
                </div>
                <div style={{display:"flex",gap:4,flexWrap:"wrap"}}>
                  {(mode==="vle"?[{id:"flash",lb:"TP Flash"},{id:"bubP",lb:"Bubble P"},{id:"bubT",lb:"Bubble T"},{id:"dewP",lb:"Dew P"},{id:"dewT",lb:"Dew T"},{id:"phflash",lb:"PH Flash"},{id:"psflash",lb:"PS Flash"}]:[{id:"txy",lb:"T-x-y"},{id:"pxy",lb:"P-x-y"}]).map(function(v){return (
                    <Pill key={v.id} active={vleType===v.id} onClick={function(){setVleType(v.id)}}>{v.lb}</Pill>
                  )})}
                </div>
              </div>
            ):null}

            {/* Conditions */}
            <div style={cardStyle}>
              <div style={{fontSize:11,color:C.textL,marginBottom:8,fontWeight:600,textTransform:"uppercase",letterSpacing:1}}>{tx("cond",lang)}</div>
              {(mode==="pure"||mode==="steam")?(
                <div style={{display:"flex",gap:4,marginBottom:10}}>
                  {(mode==="steam"?[{id:"pt",lb:"P, T"},{id:"pvf",lb:"P, VF"},{id:"tvf",lb:"T, VF"},{id:"ph",lb:"P, H"},{id:"ps",lb:"P, S"}]:[{id:"pt",lb:"P, T"},{id:"pvf",lb:"P, VF"},{id:"tvf",lb:"T, VF"}]).map(function(m){return (
                    <Pill key={m.id} active={inMode===m.id} onClick={function(){setInMode(m.id)}}>{m.lb}</Pill>
                  )})}
                </div>
              ):null}
              <div style={{display:"flex",flexDirection:"column",gap:10}}>
                {showT?(<div style={{display:"flex",alignItems:"center",gap:6}}>
                  <span style={{color:C.textL,fontSize:12,minWidth:32,fontWeight:600}}>T</span>
                  <input type="number" value={tv} onChange={function(e){setTv(e.target.value)}} style={inputStyle} />
                  <select value={tu} onChange={function(e){setTu(e.target.value)}} style={selStyle}>{TU.map(function(u){return <option key={u.id} value={u.id}>{u.lb}</option>})}</select>
                </div>):null}
                {showP?(<div style={{display:"flex",alignItems:"center",gap:6}}>
                  <span style={{color:C.textL,fontSize:12,minWidth:32,fontWeight:600}}>P</span>
                  <input type="number" value={pv} onChange={function(e){setPv(e.target.value)}} style={inputStyle} />
                  <select value={pu} onChange={function(e){setPu(e.target.value)}} style={selStyle}>{PU.map(function(u){return <option key={u.id} value={u.id}>{u.lb}</option>})}</select>
                </div>):null}
                {showVF?(<div>
                  <div style={{display:"flex",alignItems:"center",gap:6}}>
                    <span style={{color:C.textL,fontSize:12,minWidth:32,fontWeight:600}}>VF</span>
                    <input type="number" step="0.01" min="0" max="1" value={vf} onChange={function(e){setVf(e.target.value)}} placeholder="0=bub, 1=dew" style={inputStyle} />
                  </div>
                  <div style={{fontSize:10,color:C.textL,marginTop:4}}>VF=0 bubble | VF=1 dew</div>
                </div>):null}
                {showH?(<div>
                  <div style={{display:"flex",alignItems:"center",gap:6}}>
                    <span style={{color:C.textL,fontSize:12,minWidth:32,fontWeight:600}}>H</span>
                    <input type="number" value={hval} onChange={function(e){setHval(e.target.value)}} placeholder="kJ/kg or J/mol" style={inputStyle} />
                    <span style={{color:C.textL,fontSize:11}}>{mode==="vle"?"J/mol":"kJ/kg"}</span>
                  </div>
                  <div style={{fontSize:10,color:C.textL,marginTop:4}}>{mode==="vle"?(lang==="en"?"Molar enthalpy J/mol (PR EOS ref: 298.15K, 1atm)":"\u6469\u5c14\u7113 J/mol (PR\u53c2\u8003\u6001: 298.15K, 1atm)"):(lang==="en"?"Enthalpy kJ/kg (steam ~2700-3500)":"\u7113\u503c kJ/kg (\u8fc7\u70ed\u84b8\u6c7d\u7ea6 2700-3500)")}</div>
                </div>):null}
                {showS?(<div>
                  <div style={{display:"flex",alignItems:"center",gap:6}}>
                    <span style={{color:C.textL,fontSize:12,minWidth:32,fontWeight:600}}>S</span>
                    <input type="number" value={sval} onChange={function(e){setSval(e.target.value)}} placeholder="kJ/(kg·K) or J/(mol·K)" style={inputStyle} />
                    <span style={{color:C.textL,fontSize:11}}>{mode==="vle"?"J/(mol·K)":"kJ/(kg·K)"}</span>
                  </div>
                  <div style={{fontSize:10,color:C.textL,marginTop:4}}>{mode==="vle"?(lang==="en"?"Molar entropy J/(mol·K) (PR EOS ref: 298.15K, 1atm)":"\u6469\u5c14\u71b5 J/(mol\u00b7K)"):(lang==="en"?"Entropy kJ/(kg·K) (steam ~6.5-8.5)":"\u71b5\u503c kJ/(kg\u00b7K) (\u8fc7\u70ed\u84b8\u6c7d\u7ea6 6.5-8.5)")}</div>
                </div>):null}
              </div>
              <button onClick={doCalc} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#1e40af,#1d4ed8)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer",letterSpacing:0.5,boxShadow:"0 3px 12px rgba(30,64,175,0.3)",transition:"transform 0.1s"}}
                onMouseDown={function(e){e.currentTarget.style.transform="scale(0.98)"}}
                onMouseUp={function(e){e.currentTarget.style.transform="scale(1)"}}>
                {tx("calc",lang)}
              </button>
            </div>
          </div>):null}

          {/* RIGHT PANEL - Results for props / Equipment & Tools panels */}
          <div style={{flex:group==="props"?"1 1 480px":"1 1 100%",minWidth:300}}>
            {group==="props"&&err?<div style={{backgroundColor:C.errL,border:"1px solid rgba(220,38,38,0.3)",borderRadius:C.radius,padding:"10px 14px",marginBottom:12,color:C.err,fontSize:13,fontWeight:600}}>{"! "+err}</div>:null}
            {group==="props"&&res?(<div>
              {res.warn?(<div style={{padding:"6px 10px",backgroundColor:"#fff7ed",border:"1px solid #f59e0b",borderRadius:6,marginBottom:6,fontSize:11,color:"#b45309",fontWeight:600}}>{"⚠️ "}{lang==="en"?"Calculation did not fully converge (max iterations reached). Result is approximate — verify with different initial conditions.":"计算未完全收敛（达到最大迭代次数），结果为近似値——请验证或改变计算条件。"}</div>):null}
              <ResTable rows={rows} lang={lang} title={res.isVLE?"VLE | PR":(res.isSteam?"IAPWS-IF97 | H\u2082O":(meth+" | "+(res.comp?res.comp.f:"")))} />
              {res.autoNorm?(<div style={{padding:"6px 10px",backgroundColor:C.accL,border:"1px solid "+C.acc,borderRadius:6,marginTop:6,fontSize:11,color:C.acc,fontWeight:600}}>{"\u26a0\ufe0f "}{lang==="en"?"Inputs auto-normalized (sum\u22601)":"\u7ec4\u5206\u4e4b\u548c\u4e0d\u4e3a1\uff0c\u5df2\u81ea\u52a8\u5f52\u4e00\u5316 | Inputs auto-normalized (sum\u22601)"}</div>):null}
              {satC&&satC.length>0?(<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:14,marginTop:12,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:12,color:C.text,fontWeight:700,marginBottom:10}}>{tx("psatCurve",lang)} ({res.comp?res.comp.f:""})</div>
                <ResponsiveContainer width="100%" height={240}>
                  <LineChart data={satC} margin={{left:5,right:10,bottom:5,top:5}}>
                    <CartesianGrid strokeDasharray="3 3" stroke={C.border} />
                    <XAxis dataKey="T" stroke={C.borderH} tick={{fill:C.textL,fontSize:10}}
                      label={{value:"T (\u00b0C)",position:"insideBottomRight",offset:0,fill:C.textM,fontSize:10}} />
                    <YAxis stroke={C.borderH} tick={{fill:C.textL,fontSize:10}} width={45}
                      label={{value:"Psat (MPa)",angle:-90,position:"insideLeft",offset:5,fill:C.textM,fontSize:10}} />
                    <Tooltip contentStyle={{backgroundColor:C.white,border:"1px solid "+C.border,borderRadius:8,boxShadow:C.shadowM,fontSize:12}} formatter={function(v){return [typeof v==="number"?v.toFixed(4)+" MPa":v,"Psat"]}} labelFormatter={function(v){return "T = "+v+"\u00b0C"}} />
                    <Line type="monotone" dataKey="Psat" stroke={C.pri} strokeWidth={2.5} dot={false} />
                    {res.Psat?<ReferenceLine y={res.Psat/1e6} stroke={C.acc} strokeDasharray="5 5" label={{value:(res.Psat/1e6).toFixed(3),fill:C.acc,fontSize:10}} />:null}
                  </LineChart>
                </ResponsiveContainer>
              </div>):null}
            </div>):null}

            {/* Phase diagram */}
            {diagData?(<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:14,marginTop:res?12:0,boxShadow:C.shadow,border:"1px solid "+C.border}}>
              <div style={{fontSize:12,color:C.text,fontWeight:700,marginBottom:4}}>
                {diagData.type==="txy"?"T-x-y @ P="+(diagData.P/1e6).toFixed(3)+" MPa":"P-x-y @ T="+(diagData.T-273.15).toFixed(1)+"\u00b0C"}
              </div>
              <div style={{fontSize:11,color:C.textL,marginBottom:6}}>{diagData.c1.f+" (1) + "+diagData.c2.f+" (2)"}</div>
              {diagData.supercrit?(<div style={{padding:"6px 10px",backgroundColor:C.accL,border:"1px solid "+C.acc,borderRadius:6,marginBottom:8,fontSize:10,color:C.acc,lineHeight:1.5}}>
                {"\u26a0\ufe0f "}{lang==="en"
                  ?"One or more components may have limited solubility or be near-supercritical at this condition. PR EOS accuracy is reduced; some data points filtered. Verify against literature."
                  :"\u4e00\u4e2a\u6216\u591a\u4e2a\u7ec4\u5206\u5728\u6b64\u6761\u4ef6\u4e0b\u53ef\u80fd\u63a5\u8fd1\u8d85\u4e34\u754c\u6216\u6eb6\u89e3\u5ea6\u6709\u9650\uff0cPR EOS\u7cbe\u5ea6\u964d\u4f4e\u3002\u90e8\u5206\u8ba1\u7b97\u70b9\u5df2\u8fc7\u6ee4\uff0c\u8bf7\u53c2\u8003\u6587\u732e\u6838\u5b9e\u3002"}
                {diagData.Ts1!==undefined?(" Tsat("+diagData.c1.f+")="+diagData.Ts1+"\u00b0C, Tsat("+diagData.c2.f+")="+diagData.Ts2+"\u00b0C"):""}
              </div>):null}
              <div style={{display:"flex",gap:16,justifyContent:"center",marginBottom:8}}>
                <span style={{fontSize:11,color:C.pri,fontWeight:600}}>{"--- Bubble \u6ce1\u70b9 (x\u2081)"}</span>
                <span style={{fontSize:11,color:C.acc,fontWeight:600}}>{"--- Dew \u9732\u70b9 (y\u2081)"}</span>
              </div>
              <ResponsiveContainer width="100%" height={300}>
                {diagData.type==="txy"?(
                  <LineChart margin={{left:5,right:10,bottom:5,top:5}}>
                    <CartesianGrid strokeDasharray="3 3" stroke={C.border} />
                    <XAxis dataKey="x1" type="number" domain={[0,1]} stroke={C.borderH} tick={{fill:C.textL,fontSize:10}}
                      label={{value:"x\u2081, y\u2081 ("+diagData.c1.f+")",position:"insideBottomRight",offset:0,fill:C.textM,fontSize:10}} />
                    <YAxis stroke={C.borderH} tick={{fill:C.textL,fontSize:10}} width={45}
                      label={{value:"T (\u00b0C)",angle:-90,position:"insideLeft",offset:5,fill:C.textM,fontSize:10}} />
                    <Tooltip contentStyle={{backgroundColor:C.white,border:"1px solid "+C.border,borderRadius:8,boxShadow:C.shadowM,fontSize:12}} formatter={function(v){return typeof v==="number"?v.toFixed(2)+"\u00b0C":v}} />
                    <Line data={diagData.bub} type="monotone" dataKey="Tbub" stroke={C.pri} strokeWidth={2.5} dot={false} />
                    <Line data={diagData.dew.map(function(d){return{x1:d.y1,Tdew:d.Tdew}})} type="monotone" dataKey="Tdew" stroke={C.acc} strokeWidth={2.5} dot={false} />
                  </LineChart>
                ):(
                  <LineChart margin={{left:5,right:10,bottom:5,top:5}}>
                    <CartesianGrid strokeDasharray="3 3" stroke={C.border} />
                    <XAxis dataKey="x1" type="number" domain={[0,1]} stroke={C.borderH} tick={{fill:C.textL,fontSize:10}}
                      label={{value:"x\u2081, y\u2081 ("+diagData.c1.f+")",position:"insideBottomRight",offset:0,fill:C.textM,fontSize:10}} />
                    <YAxis stroke={C.borderH} tick={{fill:C.textL,fontSize:10}} width={50}
                      label={{value:"P MPa(a)",angle:-90,position:"insideLeft",offset:5,fill:C.textM,fontSize:10}} />
                    <Tooltip contentStyle={{backgroundColor:C.white,border:"1px solid "+C.border,borderRadius:8,boxShadow:C.shadowM,fontSize:12}} formatter={function(v){return typeof v==="number"?v.toFixed(4)+" MPa":v}} />
                    <Line data={diagData.bub} type="monotone" dataKey="Pbub" stroke={C.pri} strokeWidth={2.5} dot={false} />
                    <Line data={diagData.dew.map(function(d){return{x1:d.y1,Pdew:d.Pdew}})} type="monotone" dataKey="Pdew" stroke={C.acc} strokeWidth={2.5} dot={false} />
                  </LineChart>
                )}
              </ResponsiveContainer>
            </div>):null}

            {/* Empty state for props group */}
            {group==="props"&&!res&&!err&&!diagData?(<div style={{backgroundColor:C.white,borderRadius:C.radius+4,padding:48,border:"1.5px dashed "+C.border,textAlign:"center"}}>
              <div style={{fontSize:40,marginBottom:12,opacity:0.2}}>{mode==="steam"?"\ud83d\udca7":mode==="diag"?"\ud83d\udcca":"\u269b\ufe0f"}</div>
              <div style={{color:C.textL,fontSize:13}}>{mode==="steam"?tx("empty1",lang):tx("empty2",lang)}</div>
            </div>):null}

            {/* ====== EQUIPMENT: Compressor ====== */}
            {mode==="comp"?(function(){
              var iS = {flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS = {color:C.textL,fontSize:11,minWidth:70,fontWeight:600};
              var doComp = function(){
                setEqRes(null);
                var P1=parseFloat(eqIn.p1)*1e6,P2=parseFloat(eqIn.p2)*1e6,T1=parseFloat(eqIn.t1)+273.15;
                var eta=parseFloat(eqIn.eta),gam=parseFloat(eqIn.gamma),mw=parseFloat(eqIn.mw),fl=parseFloat(eqIn.flow);
                var nst=parseInt(eqIn.nst)||1;
                var tcool=parseFloat(eqIn.tcool);
                if(isNaN(P1)||isNaN(P2)||isNaN(T1)||isNaN(eta)){setErr("Check inputs");return;}
                setErr(null);
                if(nst>1){
                  if(isNaN(tcool)){setErr("Enter intercooler temperature");return;}
                  var r=calcMultiComp(P1,P2,T1,eta,gam,mw,fl>0?fl:0,nst,tcool);
                  if(r.err){setErr(r.err);return;}
                  setEqRes(Object.assign({isMulti:true},r));
                } else {
                  var r2=calcCompressor(P1,T1,P2,eta,gam,mw,fl>0?fl/3600:0);
                  if(r2.err){setErr(r2.err);return;}
                  setEqRes(r2);
                }
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:12}}>{lang==="en"?"Compressor (Isentropic)":"\u538b\u7f29\u673a (\u7b49\u71b5)"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{"T2s = T1*(P2/P1)^((\u03b3-1)/\u03b3)"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P1 MPa(a)"}</span><input type="number" value={eqIn.p1} onChange={function(e){setEq("p1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P2 MPa(a)"}</span><input type="number" value={eqIn.p2} onChange={function(e){setEq("p2",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"T1 (\u00b0C)":"T1 (\u00b0C)"}</span><input type="number" value={eqIn.t1} onChange={function(e){setEq("t1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Stages":"\u7ea7\u6570"}</span><select value={eqIn.nst} onChange={function(e){setEq("nst",e.target.value)}} style={iS}><option value="1">1</option><option value="2">2</option><option value="3">3</option><option value="4">4</option><option value="5">5</option></select></div>
                  {parseInt(eqIn.nst)>1?<div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Intercool(\u00b0C)":"\u4e2d\u51b7T(\u00b0C)"}</span><input type="number" value={eqIn.tcool} onChange={function(e){setEq("tcool",e.target.value)}} style={iS} /></div>:null}
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Eff. \u03b7s":"\u7b49\u71b5\u6548\u7387 \u03b7s"}</span><input type="number" step="0.01" value={eqIn.eta} onChange={function(e){setEq("eta",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Cp/Cv (\u03b3)"}</span><input type="number" step="0.01" value={eqIn.gamma} onChange={function(e){setEq("gamma",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>MW (g/mol)</span><input type="number" value={eqIn.mw} onChange={function(e){setEq("mw",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Mass flow(kg/h)":"\u8d28\u91cf\u6d41\u91cf(kg/h)"}</span><input type="number" value={eqIn.flow} onChange={function(e){setEq("flow",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doComp} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.isMulti?(<div style={{marginTop:12}}>
                  <ResTable rows={[{n:lang==="en"?"Stages":"\u7ea7\u6570",v:eqRes.nst,u:""},{n:lang==="en"?"Stage ratio":"\u5355\u7ea7\u538b\u6bd4",v:eqRes.rps,u:""},{n:lang==="en"?"Total Power":"\u603b\u529f\u7387",v:eqRes.W,u:"kW"}]} lang={lang} title={lang==="en"?"Multi-Stage Compressor":"\u591a\u7ea7\u538b\u7f29\u673a"} />
                  <div style={{marginTop:10,overflowX:"auto"}}>
                    <table style={{width:"100%",borderCollapse:"collapse",fontSize:11}}>
                      <thead><tr style={{backgroundColor:C.bg}}>
                        <th style={{padding:"5px 8px",textAlign:"left",color:C.textM}}>Stage</th>
                        <th style={{padding:"5px 8px",textAlign:"center",color:C.textM}}>{"Pin MPa(a)"}</th>
                        <th style={{padding:"5px 8px",textAlign:"center",color:C.textM}}>{"Pout MPa(a)"}</th>
                        <th style={{padding:"5px 8px",textAlign:"center",color:C.textM}}>{"Tin \u00b0C"}</th>
                        <th style={{padding:"5px 8px",textAlign:"center",color:C.textM}}>{"T2s \u00b0C"}</th>
                        <th style={{padding:"5px 8px",textAlign:"center",color:C.textM}}>{"T2 \u00b0C"}</th>
                        <th style={{padding:"5px 8px",textAlign:"center",color:C.textM}}>{"W kW"}</th>
                      </tr></thead>
                      <tbody>{eqRes.stages.map(function(s,i){return (<tr key={i} style={{borderBottom:"1px solid "+C.border}}>
                        <td style={{padding:"4px 8px",color:C.pri,fontWeight:700}}>{i+1}</td>
                        <td style={{padding:"4px 8px",textAlign:"center",fontFamily:"monospace"}}>{s.pin.toFixed(3)}</td>
                        <td style={{padding:"4px 8px",textAlign:"center",fontFamily:"monospace"}}>{s.pout.toFixed(3)}</td>
                        <td style={{padding:"4px 8px",textAlign:"center",fontFamily:"monospace",color:C.cTool,fontWeight:600}}>{s.tin.toFixed(1)}</td>
                        <td style={{padding:"4px 8px",textAlign:"center",fontFamily:"monospace"}}>{s.t2s.toFixed(1)}</td>
                        <td style={{padding:"4px 8px",textAlign:"center",fontFamily:"monospace",color:C.acc,fontWeight:600}}>{s.t2.toFixed(1)}</td>
                        <td style={{padding:"4px 8px",textAlign:"center",fontFamily:"monospace"}}>{s.w.toFixed(1)}</td>
                      </tr>)})}</tbody>
                    </table>
                  </div>
                </div>):null}
                {eqRes&&!eqRes.isMulti?(<div style={{marginTop:12}}>
                  <ResTable rows={[
                    {n:lang==="en"?"Compression Ratio":"\u538b\u7f29\u6bd4",v:eqRes.ratio,u:""},
                    {n:lang==="en"?"T2 isentropic":"T2 \u7b49\u71b5",v:eqRes.T2s-273.15,u:"\u00b0C"},
                    {n:lang==="en"?"T2 actual":"T2 \u5b9e\u9645",v:eqRes.T2-273.15,u:"\u00b0C"},
                    {n:lang==="en"?"Work (isentropic)":"\u7b49\u71b5\u529f",v:eqRes.ws,u:"kJ/kg"},
                    {n:lang==="en"?"Work (actual)":"\u5b9e\u9645\u529f",v:eqRes.wa,u:"kJ/kg"},
                    {n:lang==="en"?"Power":"\u529f\u7387",v:eqRes.power,u:"kW"},
                  ]} lang={lang} title={lang==="en"?"Compressor Results":"\u538b\u7f29\u673a\u8ba1\u7b97\u7ed3\u679c"} />
                </div>):null}
              </div>);
            })():null}

            {/* ====== EQUIPMENT: Expander/Turbine ====== */}
            {mode==="expd"?(function(){
              var iS = {flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS = {color:C.textL,fontSize:11,minWidth:70,fontWeight:600};
              var doExpd = function(){
                setEqRes(null);setErr(null);
                var P1=parseFloat(eqIn.p1)*1e6,P2=parseFloat(eqIn.p2)*1e6,T1=parseFloat(eqIn.t1)+273.15;
                var eta=parseFloat(eqIn.eta),gam=parseFloat(eqIn.gamma),mw=parseFloat(eqIn.mw),fl=parseFloat(eqIn.flow);
                if(isNaN(P1)||isNaN(P2)||isNaN(T1)||isNaN(eta)){setErr("Check inputs");return;}
                var r=calcExpander(P1,T1,P2,eta,gam,mw,fl>0?fl/3600:0);
                if(r.err){setErr(r.err);return;}
                setEqRes(r);
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:12}}>{lang==="en"?"Expander / Turbine (Isentropic)":"\u81a8\u80c0\u673a / \u900f\u5e73 (\u7b49\u71b5)"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:4}}>{"T2s = T1/(P1/P2)^((\u03b3-1)/\u03b3)"}</div>
                <div style={{fontSize:11,color:C.textL,marginBottom:10}}>{lang==="en"?"P1 > P2 (gas expands from high to low pressure)":"P1 > P2 (\u6c14\u4f53\u4ece\u9ad8\u538b\u81a8\u80c0\u5230\u4f4e\u538b)"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>P1 (MPa)</span><input type="number" value={eqIn.p1} onChange={function(e){setEq("p1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>P2 (MPa)</span><input type="number" value={eqIn.p2} onChange={function(e){setEq("p2",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>T1 (\u00b0C)</span><input type="number" value={eqIn.t1} onChange={function(e){setEq("t1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Eff. \u03b7s":"\u7b49\u71b5\u6548\u7387 \u03b7s"}</span><input type="number" step="0.01" value={eqIn.eta} onChange={function(e){setEq("eta",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Cp/Cv (\u03b3)</span><input type="number" step="0.01" value={eqIn.gamma} onChange={function(e){setEq("gamma",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>MW (g/mol)</span><input type="number" value={eqIn.mw} onChange={function(e){setEq("mw",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Mass flow(kg/h)":"\u8d28\u91cf\u6d41\u91cf(kg/h)"}</span><input type="number" value={eqIn.flow} onChange={function(e){setEq("flow",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doExpd} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes?(<div style={{marginTop:12}}>
                  <ResTable rows={[
                    {n:lang==="en"?"Expansion Ratio":"\u81a8\u80c0\u6bd4 P1/P2",v:eqRes.ratio,u:""},
                    {n:lang==="en"?"T2 isentropic":"T2 \u7b49\u71b5",v:eqRes.T2s-273.15,u:"\u00b0C"},
                    {n:lang==="en"?"T2 actual":"T2 \u5b9e\u9645",v:eqRes.T2-273.15,u:"\u00b0C"},
                    {n:lang==="en"?"Work output (isen.)":"\u7b49\u71b5\u8f93\u51fa\u529f",v:eqRes.ws,u:"kJ/kg"},
                    {n:lang==="en"?"Work output (actual)":"\u5b9e\u9645\u8f93\u51fa\u529f",v:eqRes.wa,u:"kJ/kg"},
                    {n:lang==="en"?"Power output":"\u8f93\u51fa\u529f\u7387",v:eqRes.power,u:"kW"},
                  ]} lang={lang} title={lang==="en"?"Expander Results":"\u81a8\u80c0\u673a\u8ba1\u7b97\u7ed3\u679c"} />
                </div>):null}
              </div>);
            })():null}

            {/* ====== EQUIPMENT: Pump ====== */}
            {mode==="pump"?(function(){
              var iS = {flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS = {color:C.textL,fontSize:11,minWidth:70,fontWeight:600};
              var doPump = function(){
                setEqRes(null);setErr(null);
                var fl=parseFloat(eqIn.flow),P1=parseFloat(eqIn.p1)*1e6,P2=parseFloat(eqIn.p2)*1e6;
                var eta=parseFloat(eqIn.eta),rho=parseFloat(eqIn.rho);
                var r=calcPump(fl,P1,P2,eta,rho);
                if(r.err){setErr(r.err);return;}
                setEqRes(r);
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:12}}>{lang==="en"?"Centrifugal Pump":"\u79bb\u5fc3\u6cf5\u8ba1\u7b97"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Flow (m\u00b3/h)":"\u6d41\u91cf (m\u00b3/h)"}</span><input type="number" value={eqIn.flow} onChange={function(e){setEq("flow",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>P1 (MPa)</span><input type="number" value={eqIn.p1} onChange={function(e){setEq("p1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>P2 (MPa)</span><input type="number" value={eqIn.p2} onChange={function(e){setEq("p2",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Eff. \u03b7":"\u6cf5\u6548\u7387 \u03b7"}</span><input type="number" step="0.01" value={eqIn.eta} onChange={function(e){setEq("eta",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Density (kg/m\u00b3)":"\u5bc6\u5ea6 (kg/m\u00b3)"}</span><input type="number" value={eqIn.rho} onChange={function(e){setEq("rho",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doPump} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes?(<div style={{marginTop:12}}>
                  <ResTable rows={[
                    {n:"\u0394P",v:eqRes.dP/1e6,u:"MPa"},
                    {n:lang==="en"?"Head":"\u626c\u7a0b",v:eqRes.head,u:"m"},
                    {n:lang==="en"?"Hydraulic Power":"\u6db2\u529b\u529f\u7387",v:eqRes.Ph,u:"kW"},
                    {n:lang==="en"?"Shaft Power":"\u8f74\u529f\u7387",v:eqRes.Ps,u:"kW"},
                  ]} lang={lang} title={lang==="en"?"Pump Results":"\u6cf5\u8ba1\u7b97\u7ed3\u679c"} />
                </div>):null}
              </div>);
            })():null}

            {/* ====== EQUIPMENT: Pipe ====== */}
            {mode==="pipe"?(function(){
              var iS = {flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS = {color:C.textL,fontSize:11,minWidth:90,fontWeight:600};
              var doPipe = function(){
                setEqRes(null);setErr(null);
                var L=parseFloat(eqIn.len),D=parseFloat(eqIn.dia),eps=parseFloat(eqIn.rough);
                var vel=parseFloat(eqIn.vel),rho=parseFloat(eqIn.rho),mu=parseFloat(eqIn.mu);
                if(isNaN(L)||isNaN(D)||isNaN(vel)||isNaN(rho)||isNaN(mu)||D<=0||mu<=0){setErr("Check inputs");return;}
                var r=calcPipe(L,D,eps,vel,rho,mu);
                setEqRes(r);
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:12}}>{lang==="en"?"Pipe Pressure Drop (Darcy-Weisbach)":"\u7ba1\u8def\u538b\u964d (Darcy-Weisbach)"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{"\u0394P = fL\u03c1v\u00b2/(2D) | Swamee-Jain f"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Length (m)":"\u7ba1\u957f L (m)"}</span><input type="number" value={eqIn.len} onChange={function(e){setEq("len",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Diameter (m)":"内径 D (m)"}</span><input type="number" value={eqIn.dia||""} onChange={function(e){var d=parseFloat(e.target.value);setEq("dia",e.target.value);var qv=parseFloat(eqIn.qv);if(qv>0&&d>0){setEq("vel",(qv/3600/(Math.PI*d*d/4)).toFixed(4));}else{var vel=parseFloat(eqIn.vel);if(vel>0&&d>0){setEq("qv",(vel*Math.PI*d*d/4*3600).toFixed(3));}}}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Roughness (m)":"\u7c97\u7cd9\u5ea6 \u03b5 (m)"}</span><input type="number" value={eqIn.rough} onChange={function(e){setEq("rough",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Vol.flow(m3/h)":"体积流量(m3/h)"}</span><input type="number" placeholder={lang==="en"?"optional":"可选"} value={eqIn.qv||""} onChange={function(e){var qv=parseFloat(e.target.value);setEq("qv",e.target.value);var d=parseFloat(eqIn.dia);if(qv>0&&d>0){setEq("vel",(qv/3600/(Math.PI*d*d/4)).toFixed(4));}}} style={Object.assign({},iS,{backgroundColor:"#f0f9ff"})} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Velocity (m/s)":"流速 (m/s)"}</span><input type="number" value={eqIn.vel||""} onChange={function(e){var vel=parseFloat(e.target.value);setEq("vel",e.target.value);var d=parseFloat(eqIn.dia);if(vel>0&&d>0){setEq("qv",(vel*Math.PI*d*d/4*3600).toFixed(3));}}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Density (kg/m\u00b3)":"\u5bc6\u5ea6 (kg/m\u00b3)"}</span><input type="number" value={eqIn.rho} onChange={function(e){setEq("rho",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Viscosity (Pa\u00b7s)":"\u7c98\u5ea6 \u03bc (Pa\u00b7s)"}</span><input type="number" value={eqIn.mu} onChange={function(e){setEq("mu",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doPipe} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes?(<div style={{marginTop:12}}>
                  <ResTable rows={[
                    {n:"Reynolds Re",v:eqRes.Re,u:""},
                    {n:lang==="en"?"Flow Regime":"\u6d41\u52a8\u72b6\u6001",v:eqRes.regime,u:""},
                    {n:lang==="en"?"Friction Factor f":"\u6469\u64e6\u56e0\u5b50 f",v:eqRes.f,u:""},
                    {n:lang==="en"?"Pressure Drop":"\u538b\u964d \u0394P",v:eqRes.dP,u:"Pa"},
                    {n:lang==="en"?"Pressure Drop":"\u538b\u964d",v:eqRes.dP/1000,u:"kPa"},
                    {n:lang==="en"?"Head Loss":"\u6c34\u5934\u635f\u5931",v:eqRes.dH,u:"m"},
                  ]} lang={lang} title={lang==="en"?"Pipe Results":"\u7ba1\u8def\u8ba1\u7b97\u7ed3\u679c"} />
                </div>):null}
                <div style={{marginTop:10,padding:8,backgroundColor:C.bg,borderRadius:6}}>
                  <div style={{fontSize:10,fontWeight:700,color:C.textM,marginBottom:4}}>{lang==="en"?"Common Fittings K-values (for Le/D or K method)":"\u5e38\u89c1\u7ba1\u4ef6\u5c40\u90e8\u963b\u529b\u7cfbK\u503c"}</div>
                  <div style={{fontSize:9,color:C.textL,fontFamily:"monospace",lineHeight:1.8}}>
                    {"90\u00b0"} {lang==="en"?"elbow (std)":"\u5f2f\u5934(\u6807\u51c6)"}: K=0.75 | {"90\u00b0"} {lang==="en"?"elbow (long)":"\u957f\u5f84\u5f2f"}: K=0.45<br/>
                    {"45\u00b0"} {lang==="en"?"elbow":"\u5f2f\u5934"}: K=0.35 | Tee({lang==="en"?"thru":"\u76f4\u901a"}): K=0.4 | Tee({lang==="en"?"branch":"\u652f\u8def"}): K=1.5<br/>
                    {lang==="en"?"Gate valve":"\u95f8\u9600"}: K=0.17 | {lang==="en"?"Globe valve":"\u622a\u6b62\u9600"}: K=6.0 | {lang==="en"?"Check valve":"\u6b62\u56de\u9600"}: K=2.0<br/>
                    {lang==="en"?"Sudden expansion":"\u7a81\u7136\u6269\u5927"}: K=1.0 | {lang==="en"?"Sudden contraction":"\u7a81\u7136\u7f29\u5c0f"}: K=0.5 | {lang==="en"?"Entrance":"\u5165\u53e3"}: K=0.5
                  </div>
                </div>
              </div>);
            })():null}
            {mode==="cop"?(function(){
              var iS = {flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS = {color:C.textL,fontSize:11,minWidth:90,fontWeight:600};
              var doCOP = function(){
                setEqRes(null);setErr(null);
                if(comps.length===0){setErr(lang==="en"?"Select a refrigerant":"\u8bf7\u9009\u62e9\u5236\u51b7\u5242");return;}
                var c0=comps[0];
                var Tevap=parseFloat(copIn.tevap)+273.15, Tcond=parseFloat(copIn.tcond)+273.15;
                var Tsh=parseFloat(copIn.tsh), Tsc=parseFloat(copIn.tsc), eta=parseFloat(copIn.eta);
                if(isNaN(Tevap)||isNaN(Tcond)||Tevap>=Tcond){setErr("Tevap < Tcond");return;}
                var r=calcRefrigCOP(Tevap,Tcond,Tsh,Tsc,eta,c0);
                if(r.err){setErr(r.err);return;}
                setEqRes(r);
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:6}}>{lang==="en"?"Refrigeration Cycle COP":"\u5236\u51b7\u5faa\u73af COP \u8ba1\u7b97"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{lang==="en"?"Select a refrigerant above, then set cycle conditions":"\u5148\u5728\u4e0a\u65b9\u9009\u62e9\u5236\u51b7\u5242\u7ec4\u5206\uff0c\u7136\u540e\u8bbe\u7f6e\u5faa\u73af\u53c2\u6570"}</div>
                <CompSearch onSelect={addComp} lang={lang} exclude={comps.map(function(c){return c.id})} />
                {comps.length>0?<div style={{marginTop:6,padding:"6px 10px",backgroundColor:C.priL,borderRadius:6,color:C.pri,fontSize:13,fontWeight:700}}>{comps[0].f+" "+comps[0].cn}</div>:null}
                <div style={{display:"flex",flexDirection:"column",gap:8,marginTop:10}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Evap T (\u00b0C)":"\u84b8\u53d1\u6e29\u5ea6 (\u00b0C)"}</span><input type="number" value={copIn.tevap} onChange={function(e){setCop("tevap",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Cond T (\u00b0C)":"\u51b7\u51dd\u6e29\u5ea6 (\u00b0C)"}</span><input type="number" value={copIn.tcond} onChange={function(e){setCop("tcond",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Superheat (\u00b0C)":"\u8fc7\u70ed\u5ea6 (\u00b0C)"}</span><input type="number" value={copIn.tsh} onChange={function(e){setCop("tsh",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Subcool (\u00b0C)":"\u8fc7\u51b7\u5ea6 (\u00b0C)"}</span><input type="number" value={copIn.tsc} onChange={function(e){setCop("tsc",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Comp. Eff. \u03b7":"\u538b\u7f29\u673a\u6548\u7387 \u03b7"}</span><input type="number" step="0.01" value={copIn.eta} onChange={function(e){setCop("eta",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doCOP} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes?(<div style={{marginTop:12}}>
                  <ResTable rows={[
                    {n:"COP",v:eqRes.COP,u:""},
                    {n:"COP (Carnot)",v:eqRes.COPcarnot,u:""},
                    {n:lang==="en"?"Carnot Eff.":"\u5361\u8bfa\u6548\u7387",v:eqRes.COP/eqRes.COPcarnot,u:""},
                    {n:"Cp/Cv (\u03b3)",v:eqRes.gamma,u:""},
                    {n:lang==="en"?"Evap P":"\u84b8\u53d1\u538b\u529b",v:eqRes.Pevap/1e6,u:"MPa"},
                    {n:lang==="en"?"Cond P":"\u51b7\u51dd\u538b\u529b",v:eqRes.Pcond/1e6,u:"MPa"},
                    {n:lang==="en"?"Press. Ratio":"\u538b\u6bd4",v:eqRes.pressRatio,u:""},
                    {n:lang==="en"?"Comp. Outlet T":"\u538b\u7f29\u673a\u51fa\u53e3T",v:eqRes.T2-273.15,u:"\u00b0C"},
                    {n:lang==="en"?"Comp. Work":"\u538b\u7f29\u529f",v:eqRes.wComp,u:"kJ/kg"},
                    {n:lang==="en"?"Evap. Capacity":"\u5236\u51b7\u91cf",v:eqRes.qEvap,u:"kJ/kg"},
                  ]} lang={lang} title={lang==="en"?"COP Results":"COP \u8ba1\u7b97\u7ed3\u679c"} />
                </div>):null}
              </div>);
            })():null}

            {/* ====== EQUIPMENT: LMTD ====== */}
            {mode==="lmtd"?(function(){
              var iS = {flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS = {color:C.textL,fontSize:11,minWidth:90,fontWeight:600};
              var doLMTD = function(){
                setEqRes(null);setErr(null);
                var r=calcLMTD(parseFloat(lmtdIn.thi),parseFloat(lmtdIn.tho),parseFloat(lmtdIn.tci),parseFloat(lmtdIn.tco),lmtdIn.flow);
                if(r.err){setErr(r.err);return;}
                setEqRes(r);
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:12}}>{lang==="en"?"LMTD (Log Mean Temperature Difference)":"LMTD \u5bf9\u6570\u5e73\u5747\u6e29\u5dee"}</div>
                <div style={{display:"flex",gap:4,marginBottom:10}}>
                  <Pill active={lmtdIn.flow==="counter"} onClick={function(){setLmtd("flow","counter")}}>{lang==="en"?"Counterflow":"\u9006\u6d41"}</Pill>
                  <Pill active={lmtdIn.flow==="parallel"} onClick={function(){setLmtd("flow","parallel")}}>{lang==="en"?"Parallel":"\u5e76\u6d41"}</Pill>
                </div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Hot In (\u00b0C)":"\u70ed\u4fa7\u8fdb\u53e3 (\u00b0C)"}</span><input type="number" value={lmtdIn.thi} onChange={function(e){setLmtd("thi",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Hot Out (\u00b0C)":"\u70ed\u4fa7\u51fa\u53e3 (\u00b0C)"}</span><input type="number" value={lmtdIn.tho} onChange={function(e){setLmtd("tho",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Cold In (\u00b0C)":"\u51b7\u4fa7\u8fdb\u53e3 (\u00b0C)"}</span><input type="number" value={lmtdIn.tci} onChange={function(e){setLmtd("tci",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Cold Out (\u00b0C)":"\u51b7\u4fa7\u51fa\u53e3 (\u00b0C)"}</span><input type="number" value={lmtdIn.tco} onChange={function(e){setLmtd("tco",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doLMTD} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes?(<div style={{marginTop:12}}>
                  <ResTable rows={[
                    {n:"LMTD",v:eqRes.LMTD,u:"\u00b0C"},
                    {n:"\u0394T1",v:eqRes.dT1,u:"\u00b0C"},
                    {n:"\u0394T2",v:eqRes.dT2,u:"\u00b0C"},
                    {n:"R (capacity ratio)",v:eqRes.R,u:""},
                    {n:"P (effectiveness)",v:eqRes.P,u:""},
                  ]} lang={lang} title="LMTD" />
                </div>):null}
              </div>);
            })():null}

            {/* ====== TOOLS: Unit Converter ====== */}
            {mode==="unit"?(function(){
              var cat = UNIT_CATS.find(function(c){return c.id===ucCat}) || UNIT_CATS[0];
              var fromU = cat.units.find(function(u){return u.id===ucFrom}) || cat.units[0];
              var toU = cat.units.find(function(u){return u.id===ucTo}) || cat.units[1];
              var val = parseFloat(ucVal);
              var result = !isNaN(val) ? toU.t(fromU.f(val)) : NaN;
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:12}}>{lang==="en"?"Unit Converter":"\u5355\u4f4d\u6362\u7b97"}</div>
                <div style={{display:"flex",gap:4,marginBottom:12,flexWrap:"wrap"}}>
                  {UNIT_CATS.map(function(c){return (
                    <Pill key={c.id} active={ucCat===c.id} onClick={function(){setUcCat(c.id);var u=c.units;setUcFrom(u[0].id);setUcTo(u.length>1?u[1].id:u[0].id)}}>
                      {c[lang]||c.en}
                    </Pill>
                  )})}
                </div>
                <div style={{display:"flex",flexDirection:"column",gap:10}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}>
                    <input type="number" value={ucVal} onChange={function(e){setUcVal(e.target.value)}} style={{flex:1,padding:"12px",backgroundColor:C.bg,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:16,fontFamily:"monospace",outline:"none"}} />
                    <select value={ucFrom} onChange={function(e){setUcFrom(e.target.value)}} style={{padding:"12px 8px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:13,cursor:"pointer"}}>
                      {cat.units.map(function(u){return <option key={u.id} value={u.id}>{u.lb}</option>})}
                    </select>
                  </div>
                  <div style={{textAlign:"center",fontSize:20,color:C.textL}}>{"\u2193"}</div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}>
                    <div style={{flex:1,padding:"12px",backgroundColor:C.priL,border:"1.5px solid "+C.pri,borderRadius:C.radius,color:C.pri,fontSize:18,fontFamily:"monospace",fontWeight:700,textAlign:"right",minHeight:24}}>
                      {!isNaN(result) ? (Math.abs(result)>1e6||Math.abs(result)<0.001&&result!==0 ? result.toExponential(6) : result.toPrecision(8)) : "—"}
                    </div>
                    <select value={ucTo} onChange={function(e){setUcTo(e.target.value)}} style={{padding:"12px 8px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:13,cursor:"pointer"}}>
                      {cat.units.map(function(u){return <option key={u.id} value={u.id}>{u.lb}</option>})}
                    </select>
                  </div>
                </div>
              </div>);
            })():null}

            {mode==="pwall"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var doPW=function(){
                setEqRes(null);setErr(null);
                var P=parseFloat(eqIn.p1);var T=parseFloat(eqIn.t1);
                if(isNaN(P)||isNaN(T)){setErr("Check P, T");return;}
                var r=calcPipeWall(P,T,pwMat,pwIdx,parseFloat(pwEw),parseFloat(pwYc),parseFloat(pwCc));
                if(r.err){setErr(r.err);return;}
                setEqRes(r);
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Pipe Wall Thickness":"\u7ba1\u9053\u58c1\u539a\u8ba1\u7b97"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:12}}>ASME B31.3 / GB/T 20801 | t = PD/(2S*E+2YP) + C</div>
                <div style={{fontSize:9,color:C.acc,marginBottom:8}}>{lang==="en"?"Note: Verify allowable stress [S] against the applicable standard before use":"\u6ce8\u610f: \u8bf7\u6838\u5b9e\u8bb8\u7528\u5e94\u529b[S]\u503c, \u4ee5\u5b9e\u9645\u6807\u51c6\u4e3a\u51c6"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Design P(MPa)":"\u8bbe\u8ba1\u538b\u529bP(MPa)"}</span><input type="number" value={eqIn.p1} onChange={function(e){setEq("p1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Design T(C)":"\u8bbe\u8ba1\u6e29\u5ea6T(C)"}</span><input type="number" value={eqIn.t1} onChange={function(e){setEq("t1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>DN</span><select value={pwIdx} onChange={function(e){setPwIdx(parseInt(e.target.value))}} style={iS}>{PIPE_SCH.map(function(p,i){return <option key={i} value={i}>{"DN"+p.dn+" (OD "+p.od+"mm)"}</option>})}</select></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Material":"\u6750\u8d28"}</span><select value={pwMat} onChange={function(e){setPwMat(e.target.value)}} style={iS}>{PIPE_MAT.map(function(m){return <option key={m.id} value={m.id}>{lang==="en"?m.en:m.cn}</option>})}</select></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>E({lang==="en"?"weld":"\u7126\u63a5"})</span><input type="number" step="0.05" value={pwEw} onChange={function(e){setPwEw(e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Y</span><input type="number" step="0.1" value={pwYc} onChange={function(e){setPwYc(e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>C({lang==="en"?"corr.mm":"\u8150\u8680mm"})</span><input type="number" step="0.5" value={pwCc} onChange={function(e){setPwCc(e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doPW} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#047857,#059669)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{lang==="en"?"Calculate":"\u8ba1\u7b97"}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.tMin?(<div style={{marginTop:12}}><ResTable rows={[{n:lang==="en"?"Standard":"\u6807\u51c6",v:eqRes.std==="ASME"?"ASME B31.3":"GB/T 20801",u:eqRes.matName},{n:"t_min",v:eqRes.tMin,u:"mm"},{n:lang==="en"?"Allowable S":"\u8bb8\u7528\u5e94\u529b[S]",v:eqRes.S,u:"MPa"},{n:"DN"+eqRes.dn,v:eqRes.od,u:"mm OD"},{n:lang==="en"?"Recommended":"\u63a8\u8350",v:eqRes.rec?eqRes.rec[0]+" ("+eqRes.rec[1]+"mm)":"N/A",u:""}].concat(eqRes.all.map(function(s){return{n:"Sch "+s[0],v:s[1],u:s[1]>=eqRes.tMin?"mm OK":"mm"}}))} lang={lang} title={(eqRes.std==="ASME"?"ASME B31.3":"GB/T 20801")} /></div>):null}
              </div>);
            })():null}

            {mode==="svsize"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var doSV=function(){
                setEqRes(null);setErr(null);
                if(svType==="gas"){var r=calcSVgas(parseFloat(svIn.w),parseFloat(svIn.t1)+273.15,parseFloat(svIn.z),parseFloat(svIn.mw),parseFloat(svIn.gamma),parseFloat(svIn.p1),parseFloat(svIn.pb),parseFloat(svIn.kd));if(r.err){setErr(r.err);return;}setEqRes(Object.assign({svt:"gas"},r))}
                else{var r2=calcSVliq(parseFloat(svIn.q),parseFloat(svIn.sg),parseFloat(svIn.p1),parseFloat(svIn.pb),parseFloat(svIn.kd));if(r2.err){setErr(r2.err);return;}setEqRes(Object.assign({svt:"liq"},r2))}
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Safety Valve Sizing":"\u5b89\u5168\u9600\u5b9a\u5f84"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>API 520 Part I | API 526 Orifice</div>
                <div style={{display:"flex",gap:4,marginBottom:10}}><Pill active={svType==="gas"} onClick={function(){setSvType("gas")}}>{lang==="en"?"Gas":"\u6c14\u4f53"}</Pill><Pill active={svType==="liq"} onClick={function(){setSvType("liq")}}>{lang==="en"?"Liquid":"\u6db2\u4f53"}</Pill></div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  {svType==="gas"?(<div style={{display:"flex",flexDirection:"column",gap:8}}>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"W (kg/h)"}</span><input type="number" value={svIn.w} onChange={function(e){setSvI("w",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"T (\u00b0C)"}</span><input type="number" value={svIn.t1} onChange={function(e){setSvI("t1",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Z</span><input type="number" value={svIn.z} onChange={function(e){setSvI("z",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"MW (g/mol)"}</span><input type="number" value={svIn.mw} onChange={function(e){setSvI("mw",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Cp/Cv (\u03b3)"}</span><input type="number" value={svIn.gamma} onChange={function(e){setSvI("gamma",e.target.value)}} style={iS} /></div>
                  </div>):(<div style={{display:"flex",flexDirection:"column",gap:8}}>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Q (L/min)"}</span><input type="number" value={svIn.q} onChange={function(e){setSvI("q",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>SG</span><input type="number" value={svIn.sg} onChange={function(e){setSvI("sg",e.target.value)}} style={iS} /></div>
                  </div>)}
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P1 kPa(a)"}</span><input type="number" value={svIn.p1} onChange={function(e){setSvI("p1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Pb kPa(a)"}</span><input type="number" value={svIn.pb} onChange={function(e){setSvI("pb",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Kd</span><input type="number" value={svIn.kd} onChange={function(e){setSvI("kd",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doSV} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{lang==="en"?"Calculate":"\u8ba1\u7b97"}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.A?(<div style={{marginTop:12}}><ResTable rows={[{n:lang==="en"?"Required A":"\u6240\u9700\u9762\u79efA",v:eqRes.A,u:"mm2"},{n:lang==="en"?"Orifice":"\u63a8\u8350\u53e3\u5f84",v:eqRes.rec?eqRes.rec[0]+" ("+eqRes.rec[1]+"mm2)":"Over",u:"API 526"}].concat(eqRes.svt==="gas"?[{n:lang==="en"?"Critical":"临界流",v:eqRes.cr?(lang==="en"?"Yes":"是"):(lang==="en"?"No":"否"),u:""},{n:"C",v:eqRes.C,u:""},{n:"Kb",v:eqRes.Kb,u:""}]:[])} lang={lang} title="API 520" /></div>):null}
              </div>);
            })():null}

            {mode==="oplate"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var doOP=function(){setEqRes(null);setErr(null);var r=calcOrifice(parseFloat(opIn.D),parseFloat(opIn.d),parseFloat(opIn.dP),parseFloat(opIn.rho),parseFloat(opIn.mu),opIn.tap);if(r.err){setErr(r.err);return;}setEqRes(r)};
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Orifice Plate Flowmeter":"\u5b54\u677f\u6d41\u91cf\u8ba1"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>ISO 5167-2:2003 | Reader-Harris/Gallagher</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>D(mm)</span><input type="number" value={opIn.D} onChange={function(e){setOp("D",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>d(mm)</span><input type="number" value={opIn.d} onChange={function(e){setOp("d",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>dP(Pa)</span><input type="number" value={opIn.dP} onChange={function(e){setOp("dP",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Density":"\u5bc6\u5ea6"}(kg/m3)</span><input type="number" value={opIn.rho} onChange={function(e){setOp("rho",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Visc.":"\u7c98\u5ea6"}(Pa.s)</span><input type="number" value={opIn.mu} onChange={function(e){setOp("mu",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Taps</span><select value={opIn.tap} onChange={function(e){setOp("tap",e.target.value)}} style={iS}><option value="flange">Flange</option><option value="corner">Corner</option><option value="D">D-D/2</option></select></div>
                </div>
                <button onClick={doOP} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#7c3aed,#8b5cf6)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{lang==="en"?"Calculate":"\u8ba1\u7b97"}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.C?(<div style={{marginTop:12}}><ResTable rows={[{n:lang==="en"?"Flow":"\u6d41\u91cf",v:eqRes.Qh,u:"m3/h"},{n:lang==="en"?"Velocity":"\u6d41\u901f",v:eqRes.v,u:"m/s"},{n:"Beta",v:eqRes.beta,u:""},{n:"C",v:eqRes.C,u:""},{n:"E",v:eqRes.E,u:""},{n:"Re",v:eqRes.Re,u:""}]} lang={lang} title="ISO 5167" /></div>):null}
              </div>);
            })():null}

            {mode==="gasflow"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:90,fontWeight:600};
              var TK=parseFloat(gfIn.t)+273.15;
              var PP=parseFloat(gfIn.p)*1e6;
              var Tstd_N=273.15, Pstd_N=101325;
              var Tstd_S=288.15, Pstd_S=101325;
              var Tstd_U=288.71, Pstd_U=101325*14.696/14.696; // 60F, 14.696psi
              // Determine which field was last edited and compute from it
              var activeField = gfIn.activeField || "am";
              var am=parseFloat(gfIn.am), nm=parseFloat(gfIn.nm), sm=parseFloat(gfIn.sm), scfm=parseFloat(gfIn.scfm), mmscfd=parseFloat(gfIn.mmscfd);
              var calcAll = function(field, val) {
                if (isNaN(TK)||isNaN(PP)||TK<=0||PP<=0) return {};
                var qa;
                if (field==="am") qa=val;
                else if (field==="nm") qa=val*(Pstd_N/PP)*(TK/Tstd_N);
                else if (field==="sm") qa=val*(Pstd_S/PP)*(TK/Tstd_S);
                else if (field==="scfm") { var us=val/35.3147*60; qa=us*(Pstd_U/PP)*(TK/Tstd_U); }
                else if (field==="mmscfd") { var ush=val*1e6/24; var us2=ush/35.3147; qa=us2*(Pstd_U/PP)*(TK/Tstd_U); }
                else return {};
                var nmv=gasStdFlow(qa,TK,PP,Tstd_N,Pstd_N);
                var smv=gasStdFlow(qa,TK,PP,Tstd_S,Pstd_S);
                var usv=gasStdFlow(qa,TK,PP,Tstd_U,Pstd_U);
                var scfhv=usv*35.3147; var scfmv=scfhv/60; var mmscfdv=scfhv*24/1e6;
                return {am:qa,nm:nmv,sm:smv,scfm:scfmv,mmscfd:mmscfdv,scfh:scfhv};
              };
              var af=gfIn.activeField||"am";
              var activeVal=parseFloat(gfIn[af]);
              var calc=(!isNaN(activeVal)&&!isNaN(TK)&&!isNaN(PP)&&TK>0&&PP>0)?calcAll(af,activeVal):{};
              var disp=function(k){if(calc[k]===undefined||isNaN(calc[k]))return"—";var v=calc[k];return (Math.abs(v)>1e6||(Math.abs(v)<0.001&&v!==0))?v.toExponential(4):v.toPrecision(6)};
              var flowFields=[
                {key:"am",lb:"Am\u00b3/h",hint:lang==="en"?"Actual flow (T, P above)":"\u5b9e\u9645\u6d41\u91cf"},
                {key:"nm",lb:"Nm\u00b3/h",hint:"0\u00b0C, 101.325 kPa"},
                {key:"sm",lb:"Sm\u00b3/h",hint:"15\u00b0C, 101.325 kPa"},
                {key:"scfm",lb:"SCFM",hint:"60\u00b0F, 14.696 psia"},
                {key:"mmscfd",lb:"MMSCFD",hint:"60\u00b0F, 14.696 psia"},
              ];
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Standard Gas Flow (Bidirectional)":"\u6807\u51c6\u6c14\u4f53\u6d41\u91cf (\u53cc\u5411\u8ba1\u7b97)"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{lang==="en"?"Q_std = Q_act \u00d7 (P/P_std) \u00d7 (T_std/T_act) | Fill any one field":"\u7406\u60f3\u6c14\u4f53\u5b9a\u5f8b | \u4efb\u610f\u586b\u5199\u4e00\u4e2a\u6d41\u91cf"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8,marginBottom:12}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"T (\u00b0C)"}</span><input type="number" value={gfIn.t} onChange={function(e){setGf(null,null,{t:e.target.value,activeField:gfIn.activeField})}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P MPa(a)"}</span><input type="number" value={gfIn.p} onChange={function(e){setGf(null,null,{p:e.target.value,activeField:gfIn.activeField})}} style={iS} /></div>
                </div>
                {(!isNaN(parseFloat(gfIn.t))&&!isNaN(parseFloat(gfIn.p))&&parseFloat(gfIn.p)>0)?null:<div style={{padding:"5px 8px",backgroundColor:C.accL,borderRadius:6,fontSize:10,color:C.acc,marginBottom:6}}>{lang==="en"?"Enter T and P(a) first":"\u8bf7\u5148\u8f93\u5165 T \u548c P(a)"}</div>}
                <div style={{display:"flex",flexDirection:"column",gap:6}}>
                  {flowFields.map(function(ff){
                    var isActive=gfIn.activeField===ff.key;
                    return (<div key={ff.key} style={{display:"flex",gap:6,alignItems:"center",padding:"6px 8px",borderRadius:8,backgroundColor:isActive?C.priL:C.bg,border:"1.5px solid "+(isActive?C.pri:C.border)}}>
                      <span style={{color:isActive?C.pri:C.textM,fontSize:11,minWidth:70,fontWeight:isActive?700:500}}>{ff.lb}</span>
                      <input type="number" value={isActive?gfIn[ff.key]:(calc[ff.key]!==undefined&&!isNaN(calc[ff.key])?parseFloat(calc[ff.key].toPrecision(6)):"")}
                        onChange={function(k){return function(e){var o={};o[k]=e.target.value;o.activeField=k;setGf(null,null,o)}}(ff.key)}
                        placeholder={isActive?"":disp(ff.key)}
                        style={{flex:1,padding:"8px 10px",backgroundColor:"transparent",border:"none",color:isActive?C.pri:C.textM,fontSize:13,fontFamily:"monospace",outline:"none"}} />
                      <span style={{fontSize:9,color:C.textL,minWidth:80,textAlign:"right"}}>{ff.hint}</span>
                    </div>);
                  })}
                </div>
                {calc.scfh!==undefined?<div style={{marginTop:8,fontSize:10,color:C.textL}}>{"SCFH = "+disp("scfh")}</div>:null}
              </div>);
            })():null}

            {mode==="mwcalc"?(function(){
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Mixture MW / Density":"\u6df7\u5408\u5206\u5b50\u91cf / \u5bc6\u5ea6"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>MW = xi*MWi | rho = PM/(RT)</div>
                <CompSearch onSelect={addComp} lang={lang} exclude={comps.map(function(c){return c.id})} />
                {comps.map(function(c){return (
                  <div key={c.id} style={{display:"flex",alignItems:"center",gap:6,marginTop:6,padding:"6px 10px",backgroundColor:C.bg,borderRadius:8}}>
                    <span style={{color:C.pri,fontWeight:700,fontFamily:"monospace",minWidth:55,fontSize:13}}>{c.f}</span>
                    <span style={{color:C.textM,fontSize:11,flex:1}}>MW={c.MW}</span>
                    <input type="number" step="0.01" placeholder="xi" value={c.mf} onChange={function(e){updMf(c.id,e.target.value)}} style={{width:58,padding:"4px 6px",backgroundColor:C.white,border:"1px solid "+C.border,borderRadius:6,color:C.text,fontSize:12,textAlign:"center",fontFamily:"monospace"}} />
                    <button onClick={function(){removeComp(c.id)}} style={{background:"none",border:"none",color:C.err,cursor:"pointer",fontSize:14,fontWeight:700}}>x</button>
                  </div>
                )})}
                {comps.length>0?(function(){
                  var fracs=comps.map(function(c){return parseFloat(c.mf)||0});
                  var sum=fracs.reduce(function(a,b){return a+b},0);
                  if(sum<=0)return <div style={{marginTop:8,color:C.textL,fontSize:11}}>{lang==="en"?"Enter mol fractions":"\u8f93\u5165\u6469\u5c14\u5206\u6570"}</div>;
                  for(var i=0;i<fracs.length;i++)fracs[i]/=sum;
                  var mwR=mixMW(comps,fracs,false);
                  var rhoS=101325*mwR.MW/(1000*R_GAS*273.15);
                  return (<div style={{marginTop:10,backgroundColor:C.bg,borderRadius:8,padding:12,fontFamily:"monospace",fontSize:13}}>
                    <div style={{display:"flex",justifyContent:"space-between",marginBottom:4}}><span style={{color:C.textM}}>Mix MW</span><span style={{color:C.pri,fontWeight:700}}>{mwR.MW.toFixed(3)} g/mol</span></div>
                    <div style={{display:"flex",justifyContent:"space-between"}}><span style={{color:C.textM}}>rho @STP</span><span style={{color:C.pri,fontWeight:700}}>{rhoS.toFixed(4)} kg/m3</span></div>
                  </div>);
                })():null}
              </div>);
            })():null}

            {mode==="cvsize"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var doCV=function(){setEqRes(null);setErr(null);if(cvType==="liq"){var r=calcCvLiq(parseFloat(cvIn.q),parseFloat(cvIn.g),parseFloat(cvIn.p1),parseFloat(cvIn.p2));if(r.err){setErr(r.err);return}setEqRes(Object.assign({cvt:"liq"},r))}else{var r2=calcCvGas(parseFloat(cvIn.w),parseFloat(cvIn.p1),parseFloat(cvIn.p2),parseFloat(cvIn.t1)+273.15,parseFloat(cvIn.mw),parseFloat(cvIn.gamma),parseFloat(cvIn.z));if(r2.err){setErr(r2.err);return}setEqRes(Object.assign({cvt:"gas"},r2))}};
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Control Valve Cv Sizing":"\u63a7\u5236\u9600Cv\u8ba1\u7b97"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>IEC 60534 / ISA S75</div>
                <div style={{display:"flex",gap:4,marginBottom:10}}><Pill active={cvType==="liq"} onClick={function(){setCvType("liq")}}>{lang==="en"?"Liquid":"\u6db2\u4f53"}</Pill><Pill active={cvType==="gas"} onClick={function(){setCvType("gas")}}>{lang==="en"?"Gas":"\u6c14\u4f53"}</Pill></div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  {cvType==="liq"?(<div style={{display:"flex",flexDirection:"column",gap:8}}>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Q(m3/h)</span><input type="number" value={cvIn.q} onChange={function(e){setCvI("q",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>SG</span><input type="number" value={cvIn.g} onChange={function(e){setCvI("g",e.target.value)}} style={iS} /></div>
                  </div>):(<div style={{display:"flex",flexDirection:"column",gap:8}}>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>W(kg/h)</span><input type="number" value={cvIn.w} onChange={function(e){setCvI("w",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>T(C)</span><input type="number" value={cvIn.t1} onChange={function(e){setCvI("t1",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>MW</span><input type="number" value={cvIn.mw} onChange={function(e){setCvI("mw",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Cp/Cv</span><input type="number" value={cvIn.gamma} onChange={function(e){setCvI("gamma",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Z</span><input type="number" value={cvIn.z} onChange={function(e){setCvI("z",e.target.value)}} style={iS} /></div>
                  </div>)}
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>P1(kPa)</span><input type="number" value={cvIn.p1} onChange={function(e){setCvI("p1",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>P2(kPa)</span><input type="number" value={cvIn.p2} onChange={function(e){setCvI("p2",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doCV} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#7c3aed,#8b5cf6)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{lang==="en"?"Calculate":"\u8ba1\u7b97"}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.Cv?(<div style={{marginTop:12}}><ResTable rows={[{n:"Cv",v:eqRes.Cv,u:""},{n:lang==="en"?"Flow Regime":"\u6d41\u52a8\u72b6\u6001",v:eqRes.regime,u:""}].concat(eqRes.cvt==="gas"?[{n:"x (dP/P1)",v:eqRes.x,u:""},{n:"Y",v:eqRes.Y,u:""},{n:"Fk",v:eqRes.Fk,u:""}]:[{n:lang==="en"?"Effective dP":"\u6709\u6548\u538b\u5dee",v:eqRes.dPeff,u:"kPa"}])} lang={lang} title="IEC 60534" /></div>):null}
              </div>);
            })():null}

            {mode==="heatval"?(function(){
              var hvIds=comps.filter(function(c){return HV_DATA[c.id]}).map(function(c){return c.id});
              var hvFracs=comps.filter(function(c){return HV_DATA[c.id]}).map(function(c){return parseFloat(c.mf)||0});
              var hvSum=hvFracs.reduce(function(a,b){return a+b},0);
              var hvR=hvSum>0?calcHeatVal(hvIds,hvFracs.map(function(f){return f/hvSum})):null;
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Fuel Gas Heating Value":"\u71c3\u6599\u6c14\u70ed\u503c\u8ba1\u7b97"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>HHV/LHV/Wobbe | ISO 6976</div>
                <CompSearch onSelect={addComp} lang={lang} exclude={comps.map(function(c){return c.id})} />
                {comps.map(function(c){return (<div key={c.id} style={{display:"flex",alignItems:"center",gap:6,marginTop:6,padding:"6px 10px",backgroundColor:C.bg,borderRadius:8}}>
                    <span style={{color:C.pri,fontWeight:700,fontFamily:"monospace",minWidth:55,fontSize:13}}>{c.f}</span>
                    <span style={{color:C.textM,fontSize:10,flex:1}}>{HV_DATA[c.id]?"HHV="+HV_DATA[c.id].HHV:"inert"}</span>
                    <input type="number" step="0.01" placeholder="xi" value={c.mf} onChange={function(e){updMf(c.id,e.target.value)}} style={{width:58,padding:"4px 6px",backgroundColor:C.white,border:"1px solid "+C.border,borderRadius:6,color:C.text,fontSize:12,textAlign:"center",fontFamily:"monospace"}} />
                    <button onClick={function(){removeComp(c.id)}} style={{background:"none",border:"none",color:C.err,cursor:"pointer",fontSize:14,fontWeight:700}}>x</button>
                  </div>)})}
                {hvR?(<div style={{marginTop:12}}><ResTable rows={[{n:"HHV",v:hvR.hhv_nm3,u:"kJ/Nm3"},{n:"LHV",v:hvR.lhv_nm3,u:"kJ/Nm3"},{n:"HHV",v:hvR.hhv_kg,u:"kJ/kg"},{n:"LHV",v:hvR.lhv_kg,u:"kJ/kg"},{n:"HHV",v:hvR.hhv_nm3/1055.06,u:"BTU/SCF"},{n:"Wobbe",v:hvR.wobbe,u:"kJ/Nm3"},{n:"MW",v:hvR.MW,u:"g/mol"},{n:lang==="en"?"Density":"\u5bc6\u5ea6",v:hvR.rho,u:"kg/Nm3"}]} lang={lang} title="ISO 6976" /></div>):<div style={{marginTop:8,color:C.textL,fontSize:11}}>{lang==="en"?"Add gas components (H2,CO,CH4...)":"\u6dfb\u52a0\u53ef\u71c3\u7ec4\u5206"}</div>}
              </div>);
            })():null}

            {mode==="vessel"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var vr=calcVessel(parseFloat(vesIn.d),parseFloat(vesIn.len),parseFloat(vesIn.h),vesIn.head,vesIn.orient);
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.pri,marginBottom:4}}>{lang==="en"?"Vessel Volume and Level":"\u5bb9\u5668\u4f53\u79ef\u4e0e\u6db2\u4f4d"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{lang==="en"?"Cylinder + head":"\u5706\u7b52 + \u5c01\u5934"}</div>
                <div style={{display:"flex",gap:4,marginBottom:10}}>
                  <Pill active={vesIn.orient==="horiz"} onClick={function(){setVesI("orient","horiz")}}>{lang==="en"?"Horizontal":"\u5367\u5f0f"}</Pill>
                  <Pill active={vesIn.orient==="vert"} onClick={function(){setVesI("orient","vert")}}>{lang==="en"?"Vertical":"\u7acb\u5f0f"}</Pill>
                </div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>D(m)</span><input type="number" value={vesIn.d} onChange={function(e){setVesI("d",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>L(m)</span><input type="number" value={vesIn.len} onChange={function(e){setVesI("len",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"h(m) from head btm":"\u6db2\u4f4dh(m) \u4ece\u5c01\u5934\u5e95"}</span><input type="number" step="0.1" value={vesIn.h} onChange={function(e){setVesI("h",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Head":"\u5c01\u5934"}</span><select value={vesIn.head} onChange={function(e){setVesI("head",e.target.value)}} style={iS}><option value="ellip">2:1 Ellipsoidal</option><option value="flat">Flat</option></select></div>
                </div>
                {vr&&!vr.err?(<div style={{marginTop:12}}><ResTable rows={[{n:lang==="en"?"Liquid Vol":"\u6db2\u4f53\u4f53\u79ef",v:vr.Vtotal,u:"m3"},{n:lang==="en"?"Liquid Vol":"\u6db2\u4f53\u4f53\u79ef",v:vr.Vtotal_L,u:"L"},{n:lang==="en"?"Total Vol":"\u603b\u5bb9\u79ef",v:vr.VtotalFull,u:"m3"},{n:lang==="en"?"Fill":"\u5145\u6ee1\u7387",v:vr.pct,u:"%"},{n:lang==="en"?"Cylinder":"\u5706\u7b52",v:vr.Vcyl,u:"m3"},{n:lang==="en"?"Heads":"\u5c01\u5934",v:vr.Vhead,u:"m3"}]} lang={lang} title={lang==="en"?"Vessel":"\u5bb9\u5668"} /></div>):null}
              </div>);
            })():null}

            {mode==="insul"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:90,fontWeight:600};
              var doIns=function(){setEqRes(null);setErr(null);var r=calcInsulation(parseFloat(insIn.d),parseFloat(insIn.tp),parseFloat(insIn.ta),parseFloat(insIn.thick),parseFloat(insIn.k),parseFloat(insIn.wind),parseFloat(insIn.emis));if(r.err){setErr(r.err);return}setEqRes(r)};
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cEquip,marginBottom:4}}>{lang==="en"?"Insulation Heat Loss":"\u4fdd\u6e29\u6563\u70ed\u8ba1\u7b97"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{lang==="en"?"Per meter of pipe":"W/m \u7ba1\u9053\u6bcf\u7c73\u6563\u70ed"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Pipe OD(mm)":"\u7ba1\u5916\u5f84(mm)"}</span><input type="number" value={insIn.d} onChange={function(e){setInsI("d",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Process T(C)":"\u5de5\u8d28\u6e29\u5ea6(C)"}</span><input type="number" value={insIn.tp} onChange={function(e){setInsI("tp",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Ambient T(C)":"\u73af\u5883\u6e29\u5ea6(C)"}</span><input type="number" value={insIn.ta} onChange={function(e){setInsI("ta",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Thick(mm, 0=bare)":"\u539a\u5ea6(mm, 0=\u88f8\u7ba1)"}</span><input type="number" value={insIn.thick} onChange={function(e){setInsI("thick",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Material":"\u6750\u6599"}</span><select value={insIn.mat} onChange={function(e){var m=INSUL_MATS.find(function(x){return x.id===e.target.value});setInsI("mat",e.target.value);if(m&&m.k>0)setInsI("k",String(m.k))}} style={iS}>{INSUL_MATS.map(function(m){return <option key={m.id} value={m.id}>{lang==="en"?m.en:m.cn}{" (k="+m.k+")"}</option>})}</select></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Conductivity k":"\u5bfc\u70ed\u7cfbk(W/mK)"}</span><input type="number" step="0.01" value={insIn.k} onChange={function(e){setInsI("k",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Wind(m/s)":"\u98ce\u901f(m/s)"}</span><input type="number" value={insIn.wind} onChange={function(e){setInsI("wind",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doIns} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,#b45309,#d97706)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{lang==="en"?"Calculate":"\u8ba1\u7b97"}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.Qbare?(<div style={{marginTop:12}}><ResTable rows={[
                  {n:lang==="en"?"Bare heat loss":"\u88f8\u7ba1\u6563\u70ed",v:eqRes.Qbare,u:"W/m"},
                  {n:lang==="en"?"Insulated loss":"\u4fdd\u6e29\u540e\u6563\u70ed",v:eqRes.Qins,u:"W/m"},
                  {n:lang==="en"?"Heat saving":"\u8282\u80fd\u7387",v:eqRes.saving,u:"%"},
                  {n:lang==="en"?"Bare surf. T":"\u88f8\u7ba1\u8868\u9762\u6e29\u5ea6",v:eqRes.Tsurf_bare,u:"C"},
                  {n:lang==="en"?"Insul. surf. T":"\u4fdd\u6e29\u8868\u9762\u6e29\u5ea6",v:eqRes.Tsurf_ins,u:"C"},
                ]} lang={lang} title={lang==="en"?"Insulation":"4fdd6e298ba17b97"} /></div>):null}
                <div style={{marginTop:8,padding:8,backgroundColor:C.bg,borderRadius:6,fontSize:9,color:C.textL,lineHeight:1.8}}>
                  {lang==="en"?"Common k values:":"5e3875285bfc70ed7cfb6570:"} Rock wool: 0.04 | Glass fiber: 0.035 | PU foam: 0.025 | Calcium silicate: 0.055 | Perlite: 0.05
                </div>
              </div>);
            })():null}

            {mode==="phasediag"?(function(){
              var phData = phaseSel==="H2O" ? WATER_PHASE : CO2_PHASE;
              var allPts = phData.getData();
              var tpT = phData.triple.T;
              var tpP = phData.triple.P;
              var cpT = phData.critical.T;
              var cpP = phData.critical.P;
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cProp,marginBottom:4}}>{"P-T"}{lang==="en"?" Phase Diagram":"\u4e09\u76f8\u56fe"}</div>
                <div style={{display:"flex",gap:4,marginBottom:10}}>
                  <Pill active={phaseSel==="H2O"} onClick={function(){setPhaseSel("H2O")}}>H2O</Pill>
                  <Pill active={phaseSel==="CO2"} onClick={function(){setPhaseSel("CO2")}}>CO2</Pill>
                </div>
                <div style={{marginBottom:6,fontSize:11,color:C.textM,lineHeight:1.8}}>
                  <span style={{color:"#059669",fontWeight:700}}>{lang==="en"?"Triple":"\u4e09\u76f8\u70b9"}: </span>{"T="+tpT+"\u00b0C, P="+tpP+" MPa"}<br/>
                  <span style={{color:"#dc2626",fontWeight:700}}>{lang==="en"?"Critical":"\u4e34\u754c\u70b9"}: </span>{"T="+cpT+"\u00b0C, P="+cpP+" MPa"}
                  {phaseSel==="CO2"?<span><br/><span style={{color:C.textL}}>{lang==="en"?"Dry ice sublimes at -78.5\u00b0C (1atm)":"\u5e72\u51b0\u5347\u534e -78.5\u00b0C (1atm)"}</span></span>:null}
                </div>
                <ResponsiveContainer width="100%" height={340}>
                  <LineChart data={allPts} margin={{top:10,right:20,bottom:25,left:15}}>
                    <CartesianGrid strokeDasharray="3 3" stroke={C.border} />
                    <XAxis dataKey="T" type="number" domain={["auto","auto"]} tick={{fontSize:10}} label={{value:"T (\u00b0C)",position:"bottom",fontSize:11,offset:0}} />
                    <YAxis type="number" scale="log" domain={["auto","auto"]} tick={{fontSize:9}} label={{value:"P MPa(a)",angle:-90,position:"insideLeft",fontSize:11,offset:5}} allowDataOverflow={true} />
                    <Tooltip contentStyle={{fontSize:10}} formatter={function(v){return typeof v==="number"?v.toFixed(4)+" MPa":v}} />
                    <Line dataKey="vap" stroke="#dc2626" dot={false} strokeWidth={2.5} connectNulls={false} name={lang==="en"?"Vaporization":"\u6c7d\u5316\u7ebf"} />
                    <Line dataKey="subl" stroke="#2563eb" dot={false} strokeWidth={2.5} connectNulls={false} name={lang==="en"?"Sublimation":"\u5347\u534e\u7ebf"} />
                    <Line dataKey="melt" stroke="#7c3aed" dot={false} strokeWidth={2.5} connectNulls={false} name={lang==="en"?"Melting":"\u7194\u5316\u7ebf"} />
                  </LineChart>
                </ResponsiveContainer>
                <div style={{display:"flex",gap:10,justifyContent:"center",fontSize:10,marginTop:6,flexWrap:"wrap"}}>
                  <span><span style={{color:"#dc2626",fontWeight:700}}>{"\u2014"}</span> {lang==="en"?"Vaporization":"\u6c7d\u5316\u7ebf(\u6db2-\u6c14)"}</span>
                  <span><span style={{color:"#2563eb",fontWeight:700}}>{"\u2014"}</span> {lang==="en"?"Sublimation":"\u5347\u534e\u7ebf(\u56fa-\u6c14)"}</span>
                  <span><span style={{color:"#7c3aed",fontWeight:700}}>{"\u2014"}</span> {lang==="en"?"Melting":"\u7194\u5316\u7ebf(\u56fa-\u6db2)"}</span>
                  {phaseSel==="H2O"?<span style={{color:C.textL,fontSize:9}}>{lang==="en"?"(Water: melting line slopes left!)":"(\u6ce8\u610f:\u6c34\u7684\u7194\u5316\u7ebf\u5411\u5de6\u503e\u659c!)"}</span>:null}
                </div>
              </div>);
            })():null}

            {mode==="ma420"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var r=null;
              if(maIn.dir==="pv2ma"){r=maConvert(parseFloat(maIn.pv),parseFloat(maIn.lo),parseFloat(maIn.hi))}
              else{r=maReverse(parseFloat(maIn.ma),parseFloat(maIn.lo),parseFloat(maIn.hi))}
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cInst,marginBottom:4}}>4-20mA {lang==="en"?"Signal Conversion":"\u4fe1\u53f7\u6362\u7b97"}</div>
                <div style={{display:"flex",gap:4,marginBottom:10}}>
                  <Pill active={maIn.dir==="pv2ma"} onClick={function(){setMaI("dir","pv2ma")}}>{lang==="en"?"PV\u2192mA":"PV\u2192mA"}</Pill>
                  <Pill active={maIn.dir==="ma2pv"} onClick={function(){setMaI("dir","ma2pv")}}>{lang==="en"?"mA\u2192PV":"mA\u2192PV"}</Pill>
                </div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Range Lo":"\u91cf\u7a0b\u4e0b\u9650"}</span><input type="number" value={maIn.lo} onChange={function(e){setMaI("lo",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Range Hi":"\u91cf\u7a0b\u4e0a\u9650"}</span><input type="number" value={maIn.hi} onChange={function(e){setMaI("hi",e.target.value)}} style={iS} /></div>
                  {maIn.dir==="pv2ma"?<div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>PV</span><input type="number" value={maIn.pv} onChange={function(e){setMaI("pv",e.target.value)}} style={iS} /></div>
                  :<div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>mA</span><input type="number" value={maIn.ma} onChange={function(e){setMaI("ma",e.target.value)}} style={iS} /></div>}
                </div>
                {r?(<div style={{marginTop:12,backgroundColor:C.bg,borderRadius:8,padding:12,fontFamily:"monospace",fontSize:13}}>
                  <div style={{display:"flex",justifyContent:"space-between",marginBottom:4}}><span style={{color:C.textM}}>{lang==="en"?"Percent":"\u767e\u5206\u6bd4"}</span><span style={{color:C.cInst,fontWeight:700}}>{r.pct.toFixed(2)}%</span></div>
                  {maIn.dir==="pv2ma"?<div style={{display:"flex",justifyContent:"space-between"}}><span style={{color:C.textM}}>mA</span><span style={{color:C.cInst,fontWeight:700}}>{r.ma.toFixed(3)} mA</span></div>
                  :<div style={{display:"flex",justifyContent:"space-between"}}><span style={{color:C.textM}}>PV</span><span style={{color:C.cInst,fontWeight:700}}>{r.pv.toFixed(4)}</span></div>}
                </div>):null}
              </div>);
            })():null}

            {mode==="tcconv"?(function(){
              var tc=TC_TYPES[tcType];
              var res=null;
              if(tc){var v=parseFloat(tcVal);if(!isNaN(v)){var isRTD=tcType==="Pt100"||tcType==="Pt1000";if(tcDir==="t2mv"){res=isRTD?{ohm:tc.toOhm(v),T:v}:{mV:tc.toMV(v),T:v}}else{res=isRTD?{ohm:v,T:tc.toT(v)}:{mV:v,T:tc.toT(v)}}}}
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cInst,marginBottom:4}}>{lang==="en"?"Thermocouple Conversion":"\u70ed\u7535\u5076\u6362\u7b97"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:8}}>ITS-90 {lang==="en"?"simplified polynomial":"\u7b80\u5316\u591a\u9879\u5f0f"}</div>
                <div style={{display:"flex",gap:4,marginBottom:8}}>
                  {["K","B","S","Pt100","Pt1000"].map(function(t){return <Pill key={t} active={tcType===t} onClick={function(){setTcType(t)}}>{TC_TYPES[t].name}</Pill>})}
                </div>
                <div style={{display:"flex",gap:4,marginBottom:10}}>
                  <Pill active={tcDir==="t2mv"} onClick={function(){setTcDir("t2mv")}}>{"T(\u00b0C)\u2192"+(tcType==="Pt100"||tcType==="Pt1000"?"\u03a9":"mV")}</Pill>
                  <Pill active={tcDir==="mv2t"} onClick={function(){setTcDir("mv2t")}}>{(tcType==="Pt100"||tcType==="Pt1000"?"\u03a9":"mV")+"\u2192T(\u00b0C)"}</Pill>
                </div>
                <div style={{display:"flex",gap:6,alignItems:"center"}}>
                  <span style={{color:C.textL,fontSize:11,minWidth:60,fontWeight:600}}>{tcDir==="t2mv"?"T(\u00b0C)":(tcType==="Pt100"||tcType==="Pt1000"?"\u03a9":"mV")}</span>
                  <input type="number" value={tcVal} onChange={function(e){setTcVal(e.target.value)}} style={{flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none"}} />
                </div>
                {res?(<div style={{marginTop:12,backgroundColor:C.bg,borderRadius:8,padding:12,fontFamily:"monospace",fontSize:13}}>
                  <div style={{display:"flex",justifyContent:"space-between",marginBottom:4}}><span style={{color:C.textM}}>{tcDir==="t2mv"?(tcType==="Pt100"||tcType==="Pt1000"?"\u03a9":"mV"):"T"}</span><span style={{color:C.cInst,fontWeight:700}}>{tcDir==="t2mv"?(res.ohm!==undefined?res.ohm.toFixed(2)+" \u03a9":res.mV.toFixed(3)+" mV"):res.T.toFixed(1)+" \u00b0C"}</span></div>
                  <div style={{display:"flex",justifyContent:"space-between"}}><span style={{color:C.textM}}>{lang==="en"?"Type":"\u7c7b\u578b"}</span><span style={{color:C.textM}}>{tc.name}</span></div>
                  <div style={{display:"flex",justifyContent:"space-between"}}><span style={{color:C.textM}}>{lang==="en"?"Range":"\u8303\u56f4"}</span><span style={{color:C.textM}}>{tc.range[0]+"\u2192"+tc.range[1]+"\u00b0C"}</span></div>
                </div>):null}
              </div>);
            })():null}

            {mode==="hxarea"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var hr=calcHXarea(parseFloat(hxIn.q),parseFloat(hxIn.u),parseFloat(hxIn.lmtd));
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cEquip,marginBottom:4}}>{lang==="en"?"Heat Exchanger Area":"\u6362\u70ed\u5668\u9762\u79ef\u4f30\u7b97"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>Q = U A {"\u0394"}Tlm</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>Q (kW)</span><input type="number" value={hxIn.q} onChange={function(e){setHxI("q",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"U (W/m\u00b2K)"}</span><input type="number" value={hxIn.u} onChange={function(e){setHxI("u",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"\u0394Tlm (\u00b0C)"}</span><input type="number" value={hxIn.lmtd} onChange={function(e){setHxI("lmtd",e.target.value)}} style={iS} /></div>
                </div>
                {hr&&!hr.err?(<div style={{marginTop:12}}><ResTable rows={[{n:lang==="en"?"Area":"\u9762\u79ef A",v:hr.A,u:"m\u00b2"},{n:"Q",v:hr.Q,u:"kW"},{n:"U",v:hr.U,u:"W/m\u00b2K"},{n:"\u0394Tlm",v:hr.LMTD,u:"\u00b0C"}]} lang={lang} title={"Q=UA\u0394T"} /></div>):null}
                <div style={{marginTop:10,fontSize:9,color:C.textL,lineHeight:1.7}}>
                  <div style={{fontWeight:700,marginBottom:2}}>{lang==="en"?"Typical U values (W/m\u00b2K)":"\u5e38\u7528U\u503c (W/m\u00b2K)"}</div>
                  {HX_U.map(function(h,i){return <div key={i}>{lang==="en"?h.hot+"/"+h.cold:h.cn}: {h.lo+"-"+h.hi+", typ "+h.typ}</div>})}
                </div>
              </div>);
            })():null}

            {mode==="npsh"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:100,fontWeight:600};
              var nr=calcNPSH(parseFloat(npshIn.ps),parseFloat(npshIn.pv),parseFloat(npshIn.hs),parseFloat(npshIn.hf),parseFloat(npshIn.rho));
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cEquip,marginBottom:4}}>NPSHa {lang==="en"?"Check":"\u6821\u6838"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>NPSHa = (Ps-Pv)/({"\u03c1"}g) + hs - hf</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Surface P (kPa)":"\u6db2\u9762\u538b\u529bPs(kPa)"}</span><input type="number" value={npshIn.ps} onChange={function(e){setNpshI("ps",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Vapor P (kPa)":"\u9971\u548c\u84b8\u6c14\u538bPv(kPa)"}</span><input type="number" value={npshIn.pv} onChange={function(e){setNpshI("pv",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Static head hs(m)":"\u9759\u538b\u5934hs(m)"}</span><input type="number" value={npshIn.hs} onChange={function(e){setNpshI("hs",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Friction hf(m)":"\u6469\u64e6\u635f\u5931hf(m)"}</span><input type="number" value={npshIn.hf} onChange={function(e){setNpshI("hf",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"\u03c1 (kg/m\u00b3)"}</span><input type="number" value={npshIn.rho} onChange={function(e){setNpshI("rho",e.target.value)}} style={iS} /></div>
                </div>
                {nr?(<div style={{marginTop:12,backgroundColor:nr.NPSHa>0?C.okL:C.errL,borderRadius:8,padding:12,textAlign:"center"}}>
                  <div style={{fontSize:11,color:C.textM}}>NPSHa</div>
                  <div style={{fontSize:28,fontWeight:700,color:nr.NPSHa>2?C.ok:C.err}}>{nr.NPSHa.toFixed(2)} m</div>
                  <div style={{fontSize:10,color:C.textM,marginTop:4}}>{nr.NPSHa>2?(lang==="en"?"OK (NPSHa > 2m)":"\u5408\u683c (NPSHa > 2m)"):(lang==="en"?"WARNING: May cavitate":"\u8b66\u544a: \u53ef\u80fd\u6c7d\u8680")}</div>
                </div>):null}
              </div>);
            })():null}

            {mode==="pipevel"?(function(){
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cTool,marginBottom:8}}>{lang==="en"?"Recommended Pipe Velocities":"\u7ba1\u9053\u63a8\u8350\u6d41\u901f"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:6}}>{lang==="en"?"Reference: Perry, GPSA, industry practice":"\u53c2\u8003: Perry/GPSA/\u884c\u4e1a\u7ecf\u9a8c"}</div>
                <div style={{overflowX:"auto"}}>
                  <table style={{width:"100%",borderCollapse:"collapse",fontSize:11}}>
                    <thead><tr style={{backgroundColor:C.bg}}>
                      <th style={{padding:"6px 8px",textAlign:"left",color:C.textM,fontWeight:600}}>{lang==="en"?"Fluid":"\u4ecb\u8d28"}</th>
                      <th style={{padding:"6px 8px",textAlign:"center",color:C.textM,fontWeight:600}}>{lang==="en"?"Low":"\u4f4e"}</th>
                      <th style={{padding:"6px 8px",textAlign:"center",color:C.cTool,fontWeight:700}}>{lang==="en"?"Typical":"\u5178\u578b"}</th>
                      <th style={{padding:"6px 8px",textAlign:"center",color:C.textM,fontWeight:600}}>{lang==="en"?"High":"\u9ad8"}</th>
                    </tr></thead>
                    <tbody>{PIPE_VEL.map(function(pv,i){return (<tr key={i} style={{borderBottom:"1px solid "+C.border}}>
                      <td style={{padding:"5px 8px",color:C.text,fontWeight:500}}>{lang==="en"?pv.fluid:pv.cn}</td>
                      <td style={{padding:"5px 8px",textAlign:"center",color:C.textM,fontFamily:"monospace"}}>{pv.lo}</td>
                      <td style={{padding:"5px 8px",textAlign:"center",color:C.cTool,fontWeight:700,fontFamily:"monospace"}}>{pv.typ}</td>
                      <td style={{padding:"5px 8px",textAlign:"center",color:C.textM,fontFamily:"monospace"}}>{pv.hi}</td>
                    </tr>)})}</tbody>
                  </table>
                </div>
                <div style={{marginTop:6,fontSize:9,color:C.textL}}>m/s</div>
              </div>);
            })():null}

            {mode==="humid"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:80,fontWeight:600};
              var hr=calcHumidity(parseFloat(humIn.t),parseFloat(humIn.rh),parseFloat(humIn.patm));
              var frostRows = hr&&hr.Tfrost!==null&&hr.Tfrost!==undefined?[{n:lang==="en"?"Frost point":"\u971c\u70b9",v:hr.Tfrost,u:"\u00b0C"}]:[];
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cTool,marginBottom:4}}>{lang==="en"?"Humidity & Dew Point":"\u6e7f\u5ea6\u4e0e\u9732\u70b9"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>Sonntag (1990) | 有效范围 -150°C 到 +100°C</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Tdb (\u00b0C)"}</span><input type="number" value={humIn.t} onChange={function(e){setHumI("t",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>RH (%)</span><input type="number" value={humIn.rh} onChange={function(e){setHumI("rh",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Patm (kPa)"}</span><input type="number" value={humIn.patm} onChange={function(e){setHumI("patm",e.target.value)}} style={iS} /></div>
                </div>
                {hr&&!isNaN(hr.Tdew)?(<div style={{marginTop:12}}><ResTable rows={[
                  {n:lang==="en"?"Dew point":"\u9732\u70b9",v:hr.Tdew,u:"\u00b0C"},
                ].concat(frostRows).concat([
                  {n:lang==="en"?"Wet bulb":"\u6e7f\u7403",v:hr.Twet,u:"\u00b0C"},
                  {n:lang==="en"?"Abs. humidity":"\u542b\u6e7f\u91cf",v:hr.W,u:"g/kg"},
                  {n:"Psat",v:hr.Psat,u:"kPa"},
                  {n:"Pw",v:hr.Pw,u:"kPa"},
                  {n:lang==="en"?"Water ppm":"\u6c34\u5206ppm",v:hr.ppm,u:"ppmv"},
                ])} lang={lang} title={lang==="en"?"Humidity":"\u6e7f\u5ea6"} /></div>):null}
                {hr&&hr.Tfrost!==null?<div style={{marginTop:6,fontSize:10,color:C.textL}}>{lang==="en"?"Frost point shown because dew point < 0\u00b0C":"\u9732\u70b9 < 0\u00b0C \u65f6\u663e\u793a\u971c\u70b9"}</div>:null}
              </div>);
            })():null}

            {/* ====== PROPS: Mixture Properties ====== */}
            {mode==="mixprop"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:100,fontWeight:600};
              var doMix=function(){
                setEqRes(null);setErr(null);
                if(!comps||comps.length<1){setErr(lang==="en"?"Select at least 1 component":"请选择至少1个组分");return;}
                var TK=parseFloat(mixIn.t)+273.15, PKPa=parseFloat(mixIn.p)*1000;
                if(isNaN(TK)||isNaN(PKPa)){setErr("Check T, P");return;}
                var z=comps.map(function(c){return parseFloat(c.frac)||0});
                var zS=z.reduce(function(a,b){return a+b},0);
                if(zS<=0){setErr(lang==="en"?"Enter mole fractions":"请输入摩尔分数");return;}
                z=z.map(function(v){return v/zS});
                var clist=comps.map(function(c){return COMPS.find(function(x){return x.id===c.id})||c});
                var res={};
                // Gas phase
                try{res.muG=gasMixViscosity(TK,z,clist)*1e6;}catch(e){}  // μPa·s
                try{res.rhoG=gasMixDensity(TK,PKPa*1000,z,clist);}catch(e){}
                // Liquid phase
                try{var lv=liqMixViscosity(TK,z,clist);res.muL=lv.mu*1000;res.warnLLE=lv.warn;}catch(e){}
                try{res.rhoL=liqMixDensity(TK,z,clist);}catch(e){}
                // MW mix
                var MWmix=0;for(var i=0;i<clist.length;i++)MWmix+=z[i]*clist[i].MW;
                res.MW=MWmix; res.T=TK; res.P=PKPa;
                setEqRes(res);
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cProp,marginBottom:4}}>{lang==="en"?"Mixture Properties":"混合物性"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{lang==="en"?"Gas: Lucas+Herning-Zipperer | Liquid: Letsou-Stiel+Lobe | Density: DIPPR-105+Rackett":"气: Lucas+Herning-Zipperer | 液: Letsou-Stiel+Lobe | 密度: DIPPR-105+Rackett"}</div>
                <div style={{fontSize:11,color:C.textM,marginBottom:8,fontWeight:600}}>{lang==="en"?"Use VLE component selector above to set composition":"使用上方VLE组分选择器设置组成"}</div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"T (°C)"}</span><input type="number" value={mixIn.t} onChange={function(e){setMixI("t",e.target.value)}} style={iS} /></div>
                  <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P (MPa)"}</span><input type="number" value={mixIn.p} onChange={function(e){setMixI("p",e.target.value)}} style={iS} /></div>
                </div>
                <button onClick={doMix} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,"+C.cProp+",#7c3aed)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.MW?(
                  <div style={{marginTop:12}}>
                    {eqRes.warnLLE?<div style={{padding:"6px 10px",backgroundColor:"#fff7ed",border:"1px solid #f59e0b",borderRadius:6,marginBottom:8,fontSize:10,color:"#b45309"}}>
                      {"⚠️ "}{lang==="en"?"H₂O + hydrocarbon system detected. Liquid phase may split — Lobe mixing rule may underestimate viscosity.":"检测到H₂O+烃类体系，液相可能分层，Lobe混合规则粘度可能偏低。"}
                    </div>:null}
                    <ResTable rows={[
                      {n:lang==="en"?"MW mix":"混合MW",v:eqRes.MW,u:"g/mol"},
                      {n:lang==="en"?"Gas viscosity μ":"气相粘度 μ",v:eqRes.muG,u:"μPa·s"},
                      {n:lang==="en"?"Gas density ρ":"气相密度 ρ",v:eqRes.rhoG,u:"kg/m³"},
                      {n:lang==="en"?"Liq viscosity μ":"液相粘度 μ",v:eqRes.muL,u:"mPa·s"},
                      {n:lang==="en"?"Liq density ρ":"液相密度 ρ",v:eqRes.rhoL,u:"kg/m³"},
                    ]} lang={lang} title={lang==="en"?"Mixture Properties":"混合物物性"} />
                    <div style={{marginTop:6,fontSize:9,color:C.textL}}>{lang==="en"?"Liquid values valid below critical T. Gas viscosity: Lucas method. Density: DIPPR-105/Rackett ideal mixing.":"液相值在临界温度以下有效。气相粘度: Lucas法。密度: DIPPR-105/Rackett理想混合。"}</div>
                  </div>
                ):null}
              </div>);
            })():null}

            {/* ====== INST: Burst Disc (爆破片) ====== */}
            {mode==="burstdisc"?(function(){
              var iS={flex:1,padding:"10px 12px",backgroundColor:C.white,border:"1.5px solid "+C.border,borderRadius:C.radius,color:C.text,fontSize:14,fontFamily:"monospace",outline:"none",boxSizing:"border-box"};
              var lS={color:C.textL,fontSize:11,minWidth:110,fontWeight:600};
              var bdType=bdIn.type, bdStd=bdIn.std;
              var doBD=function(){
                setEqRes(null);setErr(null);
                if(bdType==="gas"){
                  var W=parseFloat(bdIn.w),T=parseFloat(bdIn.t1)+273.15,Z=parseFloat(bdIn.z),MW=parseFloat(bdIn.mw),g=parseFloat(bdIn.gamma),P0=parseFloat(bdIn.p0),Kd=parseFloat(bdIn.kd);
                  if([W,T,Z,MW,g,P0,Kd].some(isNaN)||P0<=0){setErr("Check inputs");return;}
                  var r=calcBurstDiscGas(W,T,Z,MW,g,P0,bdStd,Kd);
                  setEqRes(Object.assign({bdt:"gas"},r));
                } else {
                  var Q=parseFloat(bdIn.q),G=parseFloat(bdIn.sg),P1=parseFloat(bdIn.p1),P2=parseFloat(bdIn.p2),Kd2=parseFloat(bdIn.kd),Kw=parseFloat(bdIn.kw),Kv=parseFloat(bdIn.kv);
                  if([Q,G,P1,P2].some(isNaN)){setErr("Check inputs");return;}
                  var r2=calcBurstDiscLiq(Q,G,P1,P2,bdStd,Kd2,Kw,Kv);
                  if(r2.err){setErr(r2.err);return;}
                  setEqRes(Object.assign({bdt:"liq"},r2));
                }
              };
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cInst,marginBottom:4}}>{lang==="en"?"Burst Disc Sizing":"爆破片面积计算"}</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:10}}>{lang==="en"?"GB/T 567.2-2012 | ISO 4126-2:2003":"GB/T 567.2-2012 | ISO 4126-2:2003"}</div>
                {/* Standard selector */}
                <div style={{display:"flex",gap:8,marginBottom:10}}>
                  {["GB","ISO"].map(function(s){return <button key={s} onClick={function(){setBd("std",s)}} style={{flex:1,padding:"8px",backgroundColor:bdStd===s?C.cInst:"transparent",border:"1.5px solid "+(bdStd===s?C.cInst:C.border),borderRadius:C.radius,color:bdStd===s?"#fff":C.textM,fontWeight:700,cursor:"pointer",fontSize:12}}>{s==="GB"?"GB/T 567.2":"ISO 4126-2"}</button>})}
                </div>
                {/* Type selector */}
                <div style={{display:"flex",gap:8,marginBottom:12}}>
                  {["gas","liq"].map(function(t){return <button key={t} onClick={function(){setBd("type",t)}} style={{flex:1,padding:"8px",backgroundColor:bdType===t?C.cInst:"transparent",border:"1.5px solid "+(bdType===t?C.cInst:C.border),borderRadius:C.radius,color:bdType===t?"#fff":C.textM,fontWeight:600,cursor:"pointer",fontSize:12}}>{t==="gas"?(lang==="en"?"Gas/Vapor":"气体/蒸汽"):(lang==="en"?"Liquid":"液体")}</button>})}
                </div>
                <div style={{display:"flex",flexDirection:"column",gap:8}}>
                  {bdType==="gas"?(<>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"W (kg/h)":"质量流量 W (kg/h)"}</span><input type="number" value={bdIn.w} onChange={function(e){setBd("w",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"T1 (°C)":"温度 T1 (°C)"}</span><input type="number" value={bdIn.t1} onChange={function(e){setBd("t1",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"Z (compress.)":"压缩因子 Z"}</span><input type="number" value={bdIn.z} onChange={function(e){setBd("z",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"MW (g/mol)":"分子量 MW"}</span><input type="number" value={bdIn.mw} onChange={function(e){setBd("mw",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"γ (Cp/Cv)":"比热容 γ"}</span><input type="number" value={bdIn.gamma} onChange={function(e){setBd("gamma",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P0 kPa(a)"}</span><input type="number" value={bdIn.p0} onChange={function(e){setBd("p0",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Kd (默认0.62)"}</span><input type="number" value={bdIn.kd} onChange={function(e){setBd("kd",e.target.value)}} style={iS} /></div>
                  </>):(<>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Q (L/min)"}</span><input type="number" value={bdIn.q} onChange={function(e){setBd("q",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{lang==="en"?"SG (密度比)":"密度比 G"}</span><input type="number" value={bdIn.sg} onChange={function(e){setBd("sg",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P1 kPa(a)"}</span><input type="number" value={bdIn.p1} onChange={function(e){setBd("p1",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"P2 kPa(a)"}</span><input type="number" value={bdIn.p2} onChange={function(e){setBd("p2",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Kd (默认0.62)"}</span><input type="number" value={bdIn.kd} onChange={function(e){setBd("kd",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Kw (默认1.0)"}</span><input type="number" value={bdIn.kw} onChange={function(e){setBd("kw",e.target.value)}} style={iS} /></div>
                    <div style={{display:"flex",gap:6,alignItems:"center"}}><span style={lS}>{"Kv (默认1.0)"}</span><input type="number" value={bdIn.kv} onChange={function(e){setBd("kv",e.target.value)}} style={iS} /></div>
                  </>)}
                </div>
                <button onClick={doBD} style={{width:"100%",marginTop:14,padding:"12px",background:"linear-gradient(135deg,"+C.cInst+",#0369a1)",border:"none",borderRadius:C.radius,color:"#fff",fontSize:15,fontWeight:700,cursor:"pointer"}}>{tx("calc",lang)}</button>
                {err?<div style={{marginTop:8,color:C.err,fontSize:12}}>{"! "+err}</div>:null}
                {eqRes&&eqRes.A?(
                  <div style={{marginTop:12}}>
                    <ResTable rows={[
                      {n:lang==="en"?"Required Area A":"所需面积 A",v:eqRes.A,u:"mm²"},
                      {n:"Kd",v:eqRes.Kd,u:""},
                      eqRes.bdt==="gas"?{n:"C (气体流量系数)",v:eqRes.C,u:""}:{n:"实际压差 ΔP",v:eqRes.dP,u:"kPa"},
                    ].filter(Boolean)} lang={lang} title={bdStd+" 爆破片"} />
                    <div style={{marginTop:8,padding:"8px 10px",backgroundColor:"#f0f9ff",border:"1px solid #bae6fd",borderRadius:6,fontSize:10,color:"#0369a1",lineHeight:1.6}}>
                      {"⚠️ "}{lang==="en"?"This is a preliminary calculated area. Final area must be verified by the bursting disc manufacturer using their certified Kd value per the applicable standard (ASME UV / 国内TS认证).":"以上为初步计算面积。最终面积必须由供应商按其认证Kd值根据相应标准（ASME UV/国内TS认证）进行反算核验确认。"}
                    </div>
                  </div>
                ):null}
              </div>);
            })():null}


            {mode==="b3610"?(function(){
              return (<div style={{backgroundColor:C.white,borderRadius:C.radius,padding:16,boxShadow:C.shadow,border:"1px solid "+C.border}}>
                <div style={{fontSize:14,fontWeight:700,color:C.cTool,marginBottom:4}}>ASME B36.10M</div>
                <div style={{fontSize:10,color:C.textL,marginBottom:8}}>{lang==="en"?"Carbon & Alloy Steel Pipe":"\u78b3\u94a2/\u5408\u91d1\u94a2\u7ba1"}</div>
                <div style={{overflowX:"auto"}}>
                  <table style={{width:"100%",borderCollapse:"collapse",fontSize:10}}>
                    <thead><tr style={{backgroundColor:C.bg}}>
                      <th style={{padding:"4px 6px",textAlign:"left",color:C.textM}}>NPS</th>
                      <th style={{padding:"4px 6px",textAlign:"center",color:C.textM}}>DN</th>
                      <th style={{padding:"4px 6px",textAlign:"center",color:C.textM}}>OD</th>
                      <th style={{padding:"4px 6px",textAlign:"center",color:C.textM}}>STD</th>
                      <th style={{padding:"4px 6px",textAlign:"center",color:C.textM}}>XS</th>
                      <th style={{padding:"4px 6px",textAlign:"center",color:C.textM}}>160</th>
                      <th style={{padding:"4px 6px",textAlign:"center",color:C.textM}}>XXS</th>
                    </tr></thead>
                    <tbody>{B3610.map(function(r,i){return (<tr key={i} style={{borderBottom:"1px solid "+C.border}}>
                      <td style={{padding:"3px 6px",fontWeight:600,color:C.text}}>{r.nps}</td>
                      <td style={{padding:"3px 6px",textAlign:"center",color:C.textM}}>{r.dn}</td>
                      <td style={{padding:"3px 6px",textAlign:"center",color:C.pri,fontFamily:"monospace"}}>{r.od}</td>
                      <td style={{padding:"3px 6px",textAlign:"center",fontFamily:"monospace"}}>{r.sch.STD||"-"}</td>
                      <td style={{padding:"3px 6px",textAlign:"center",fontFamily:"monospace"}}>{r.sch.XS||"-"}</td>
                      <td style={{padding:"3px 6px",textAlign:"center",fontFamily:"monospace"}}>{r.sch["160"]||"-"}</td>
                      <td style={{padding:"3px 6px",textAlign:"center",fontFamily:"monospace"}}>{r.sch.XXS||"-"}</td>
                    </tr>)})}</tbody>
                  </table>
                </div>
                <div style={{marginTop:4,fontSize:8,color:C.textL}}>{lang==="en"?"Wall thickness (mm)":"\u58c1\u539a(mm)"}</div>
              </div>);
            })():null}

          </div>
        </div>

        {showInfo?(<div onClick={function(){setShowInfo(false)}} style={{position:"fixed",top:0,left:0,right:0,bottom:0,backgroundColor:"rgba(0,0,0,0.5)",zIndex:9999,display:"flex",alignItems:"center",justifyContent:"center",padding:20}}>
          <div onClick={function(e){e.stopPropagation()}} style={{backgroundColor:C.white,borderRadius:12,padding:24,maxWidth:540,width:"100%",maxHeight:"80vh",overflowY:"auto",boxShadow:"0 20px 60px rgba(0,0,0,0.3)"}}>
            <div style={{fontSize:20,fontWeight:700,color:C.pri,marginBottom:4}}>ChemCalc</div>
            <div style={{fontSize:12,color:C.textM,marginBottom:16}}>{lang==="en"?"Chemical Engineering Calculator":"\u5316\u5b66\u5de5\u7a0b\u8ba1\u7b97\u5e73\u53f0"}</div>
            <div style={{fontSize:11,color:C.text,lineHeight:1.8}}>
              <div style={{fontWeight:700,color:C.cProp,marginTop:4,fontSize:12}}>{lang==="en"?"PROPERTIES (6 modules)":"\u7269\u6027\u8ba1\u7b97 (6\u4e2a\u6a21\u5757)"}</div>
              <div style={{color:C.textM,fontSize:10}}>
                <b>Pure:</b>{" "}{lang==="en"?"Z, ρ, H, S, Cp, Cv, Psat, Vm (PR/SRK/Ideal), DIPPR-105 liquid density":"Z, \u5bc6\u5ea6, H, S, Cp, Cv, Psat, Vm (PR/SRK/\u7406\u60f3\u6c14\u4f53), DIPPR-105\u6db2\u76f8\u5bc6\u5ea6"}<br/>
                <b>Steam:</b>{" "}{lang==="en"?"IAPWS-IF97 R1/R2/R4: P-T, P-VF, T-VF, P-H, P-S flash (Region 3 boundary check)":"IAPWS-IF97 \u533a\u57df1/2/4\uff1aP-T, P-VF, T-VF, P-H, P-S\u95ea\u84b8 (\u533a\u57df3\u8fb9\u754c\u68c0\u67e5)"}<br/>
                <b>VLE:</b>{" "}{lang==="en"?"TP/PH/PS Flash, Bubble P/T, Dew P/T (PR EOS + successive substitution); mixture H, S":"TP/PH/PS\u95ea\u84b8, \u6ce1\u70b9P/T, \u9732\u70b9P/T (PR EOS)\uff1b\u6df7\u5408\u7269H, S"}<br/>
                <b>Diagram:</b>{" T-x-y, P-x-y (PR EOS)"}<br/>
                <b>{lang==="en"?"Mix Props":"\u6df7\u5408\u7269\u6027"}</b>{" "}{lang==="en"?"μ gas (Lucas+H-Z), μ liq (Letsou-Stiel+Lobe), ρ gas/liq (PR EOS/DIPPR-105)":"\u6c14\u76f8\u7c98\u5ea6 Lucas+H-Z\uff1b\u6db2\u76f8\u7c98\u5ea6 Letsou-Stiel+Lobe\uff1b\u6c14/\u6db2\u76f8\u5bc6\u5ea6"}
              </div>
              <div style={{color:C.textL,fontSize:9,marginTop:2}}>{"36 \u7ec4\u5206: H\u2082/N\u2082/O\u2082/Ar/He/Ne | CO/CO\u2082/H\u2082O/H\u2082S/COS | CH\u2084-nC\u2084/\u82ef/\u7532\u82ef/\u4e8c\u7532\u82ef/\u73af\u5df1\u70f7 | MeOH/EtOH/\u4e19\u916e/NH\u2083/SO\u2082/HCl/Cl\u2082 | R22/R32/R134a/R125/R152a/R1234yf | 70+ k\u1d62\u2c7c"}</div>

              <div style={{fontWeight:700,color:C.cEquip,marginTop:8,fontSize:12}}>{lang==="en"?"EQUIPMENT (8 modules)":"\u8bbe\u5907\u8ba1\u7b97 (8\u4e2a\u6a21\u5757)"}</div>
              <div style={{color:C.textM,fontSize:10}}>
                <b>{lang==="en"?"Compressor":"\u538b\u7f29\u673a"}</b>{" "}{lang==="en"?"Isentropic, 1-5 stages, intercooling, per-stage T/W":"\u7b49\u71b5, 1-5\u7ea7, \u4e2d\u95f4\u51b7\u5374, \u6bcf\u7ea7T/W"}<br/>
                <b>{lang==="en"?"Expander":"\u81a8\u80c0\u673a"}</b>{" "}{lang==="en"?"Isentropic turbine/expander":"\u7b49\u71b5\u900f\u5e73/\u81a8\u80c0\u673a"}<br/>
                <b>Pump</b>{" "}{lang==="en"?"Hydraulic power, head | NPSH check":"\u6db2\u529b\u529f\u7387, \u626c\u7a0b | NPSH\u6838\u9a8c"}<br/>
                <b>{lang==="en"?"Pipe \u0394P":"\u7ba1\u9053\u538b\u964d"}</b>{" Darcy-Weisbach + Swamee-Jain"}<br/>
                <b>{lang==="en"?"Refrig. COP":"\u5236\u51b7COP"}</b>{" "}{lang==="en"?"Vapor compression cycle, 5 refrigerants":"\u4e94\u79cd\u5236\u51b7\u5242\u84b8\u6c14\u538b\u7f29\u5faa\u73af"}<br/>
                <b>LMTD</b>{" "}{lang==="en"?"Counter/parallel + HX area":"\u9006\u6d41/\u5e76\u6d41 + \u6362\u70ed\u9762\u79ef"}<br/>
                <b>{lang==="en"?"Insulation":"\u4fdd\u6e29\u6563\u70ed"}</b>{" "}{lang==="en"?"Pipe/vessel heat loss, 7 insulation materials":"\u7ba1\u9053/\u5bb9\u5668\u6563\u70ed, 7\u79cd\u4fdd\u6e29\u6750\u6599"}<br/>
                <b>HX Area</b>{" Q = U·A·LMTD, 10 fluid pairs"}
              </div>

              <div style={{fontWeight:700,color:C.cTool,marginTop:8,fontSize:12}}>{lang==="en"?"TOOLS (9 modules)":"\u5de5\u7a0b\u5de5\u5177 (9\u4e2a\u6a21\u5757)"}</div>
              <div style={{color:C.textM,fontSize:10}}>
                <b>{lang==="en"?"Unit Conv.":"\u5355\u4f4d\u6362\u7b97"}</b>{" 10 "}{lang==="en"?"categories, auto-convert on unit switch":"\u7c7b, \u5207\u6362\u5355\u4f4d\u81ea\u52a8\u8f6c\u6362"}<br/>
                <b>{lang==="en"?"Pipe Wall":"\u7ba1\u9053\u58c1\u539a"}</b>{" ASME B31.3 / GB/T 20801, 10 "}{lang==="en"?"materials":"种材料"}<br/>
                <b>{lang==="en"?"Std Flow":"\u6807\u51c6\u6d41\u91cf"}</b>{" "}{lang==="en"?"Bidirectional: Am³/h ↔ Nm³/h ↔ Sm³/h ↔ SCFM ↔ MMSCFD":"\u53cc\u5411: Am\u00b3/h \u2194 Nm\u00b3/h \u2194 Sm\u00b3/h \u2194 SCFM \u2194 MMSCFD"}<br/>
                <b>{lang==="en"?"Orifice":"\u5b54\u677f"}</b>{" ISO 5167-2, Reader-Harris/Gallagher"}<br/>
                <b>{lang==="en"?"Heat Value":"\u70ed\u503c"}</b>{" HHV/LHV/Wobbe (ISO 6976), 12 "}{lang==="en"?"components":"组分"}<br/>
                <b>{lang==="en"?"Vessel Vol.":"\u5bb9\u5668\u6db2\u4f4d"}</b>{" "}{lang==="en"?"H/V cylinder + 2:1 ellipsoidal heads":"\u6c34\u5e73/\u5782\u76f4\u5706\u7b52 + 2:1\u692d\u5706\u5c01\u5934"}<br/>
                <b>MW Calc</b>{" "}{lang==="en"?"Mixture MW, density, mol/mass fraction":"\u6df7\u5408\u5206\u5b50\u91cf, \u5bc6\u5ea6, \u6469\u5c14/\u8d28\u91cf\u5206\u6570"}<br/>
                <b>{lang==="en"?"Pipe Sch.":"\u7ba1\u89c4\u8868"}</b>{" ASME B36.10M"}<br/>
                <b>{lang==="en"?"Velocity Ref.":"\u6d41\u901f\u53c2\u8003"}</b>{" 14 "}{lang==="en"?"fluid types":"流体类型"}
              </div>

              <div style={{fontWeight:700,color:C.cInst,marginTop:8,fontSize:12}}>{lang==="en"?"INSTRUMENTS (6 modules)":"\u4eea\u8868\u81ea\u63a7 (6\u4e2a\u6a21\u5757)"}</div>
              <div style={{color:C.textM,fontSize:10}}>
                <b>4-20mA</b>{" "}{lang==="en"?"PV↔signal conversion, custom range":"\u6d4b\u91cf\u4fe1\u53f7\u4e92\u8f6c, \u81ea\u5b9a\u4e49\u91cf\u7a0b"}<br/>
                <b>T/C & RTD</b>{" Pt100/Pt1000/K/B/S, T↔mV/\u03a9"}<br/>
                <b>{lang==="en"?"Orifice":"\u5b54\u677f"}</b>{" ISO 5167-2 (also in Tools)"}<br/>
                <b>{lang==="en"?"Valve Cv":"\u9600\u95e8Cv"}</b>{" IEC 60534, liquid + gas (N8=34.7 corrected)"}<br/>
                <b>{lang==="en"?"Safety Valve":"\u5b89\u5168\u9600"}</b>{" API 520/526, gas + liquid"}<br/>
                <b>{lang==="en"?"Burst Disc":"\u7206\u7834\u7247"}</b>{" GB/T 567.2-2012 + ISO 4126-2, gas + liquid"}
              </div>

              <div style={{marginTop:10,padding:8,backgroundColor:C.bg,borderRadius:6}}>
                <div style={{fontSize:10,fontWeight:700,color:C.textM,marginBottom:3}}>{lang==="en"?"v0.1 Changes":"\u672c\u7248\u4e3b\u8981\u4fee\u590d"}</div>
                <div style={{fontSize:9,color:C.textL,lineHeight:1.8}}>
                  {"• "}{lang==="en"?"Liquid Cv: corrected formula Q/N1×√(G/ΔP), matches IEC 60534":"液体Cv公式修正，与IEC 60534完全一致"}<br/>
                  {"• "}{lang==="en"?"Gas Cv: N8 corrected from 94.8 to 34.7, results now match industry tools":"气体Cv N8从94.8修正为34.7，结果与行业工具吻合"}<br/>
                  {"• "}{lang==="en"?"New: Mixture properties (μ, ρ gas/liquid, Lucas+Herning-Zipperer+Lobe)":"新增：混合物物性计算（气液粘度、密度）"}<br/>
                  {"• "}{lang==="en"?"New: Burst disc sizing (GB/T 567.2-2012 + ISO 4126-2)":"新增：爆破片面积计算（GB/T 567.2 + ISO 4126-2）"}<br/>
                  {"• "}{lang==="en"?"New: B-type and S-type thermocouples (ITS-90, NIST verified)":"新增：B型和S型热电偶（ITS-90，经NIST数据验证）"}<br/>
                  {"• "}{lang==="en"?"Safety valve moved to Instruments category":"安全阀归入仪表自控分类"}<br/>
                  {"• "}{lang==="en"?"Pipe: vol.flow ↔ velocity bidirectional linkage":"管道压降：体积流量与流速双向联动"}<br/>
                  {"• "}{lang==="en"?"Dew point range extended to -150°C":"露点计算范围扩展至-150°C"}<br/>
                  {"• "}{lang==="en"?"Humidity formula display fixed (garbled unicode)":"湿度公式显示乱码修复"}
                </div>
              </div>

              <div style={{marginTop:8,padding:8,backgroundColor:C.bg,borderRadius:6}}>
                <div style={{fontSize:10,fontWeight:700,color:C.textM,marginBottom:2}}>{lang==="en"?"Standards & References":"\u6807\u51c6\u4e0e\u6570\u636e\u6765\u6e90"}</div>
                <div style={{fontSize:9,color:C.textL,lineHeight:1.7}}>{"IAPWS-IF97 (2007) | Peng-Robinson (1976) | SRK (1972) | DIPPR 105 | ASME B31.3 | GB/T 20801 | API 520/526 | ISO 5167-2:2003 | IEC 60534 | ISO 6976 | Sonntag (1990) | ASME B36.10M | Perry's 8th Ed."}</div>
              </div>

              <div style={{marginTop:8,padding:8,backgroundColor:C.errL,borderRadius:6,border:"1px solid rgba(220,38,38,0.15)"}}>
                <div style={{fontSize:9,color:C.textM,lineHeight:1.6}}>
                  {lang==="en"
                    ?"Disclaimer: Data verified where possible; errors may exist. Always verify independently before engineering use. Developer assumes no liability."
                    :"\u58f0\u660e\uff1a\u672c\u5de5\u5177\u5df2\u5c3d\u529b\u6838\u5b9e\u6570\u636e\uff0c\u4f46\u4ecd\u53ef\u80fd\u5b58\u5728\u8bef\u5dee\u3002\u5b9e\u9645\u5de5\u7a0b\u5e94\u7528\u524d\u52a1\u5fc5\u518d\u6b21\u6838\u5b9e\uff0c\u5f00\u53d1\u8005\u4e0d\u627f\u62c5\u4f7f\u7528\u8d23\u4efb\u3002"}
                </div>
              </div>
              <div style={{marginTop:6,fontSize:9,color:C.textL,textAlign:"center"}}>
                {"\u6b22\u8fce\u53cd\u9988 / Feedback: chemcalc@outlook.com"}
              </div>
            </div>
            <button onClick={function(){setShowInfo(false)}} style={{width:"100%",marginTop:16,padding:"10px",background:C.pri,border:"none",borderRadius:8,color:"#fff",fontSize:14,fontWeight:600,cursor:"pointer"}}>{lang==="en"?"Close":"\u5173\u95ed"}</button>
          </div>
        </div>):null}

        {/* Footer */}
        <div style={{marginTop:24,padding:"16px 0",borderTop:"1px solid "+C.border,textAlign:"center"}}>
          <div style={{color:C.textL,fontSize:10,marginBottom:4}}>
            ChemCalc v0.1 | 36 Components | PR/SRK/IAPWS-IF97 | ASME B31.3 | API 520 | ISO 5167
          </div>
          <div style={{color:C.textL,fontSize:9,marginBottom:4,lineHeight:1.5}}>
            {lang==="en"
              ?"Disclaimer: Data has been verified where possible, but errors may exist. Always verify independently before use in real engineering. The developer assumes no liability."
              :"\u58f0\u660e\uff1a\u672c\u5de5\u5177\u5df2\u5c3d\u529b\u6838\u5b9e\u6570\u636e\uff0c\u4f46\u4ecd\u53ef\u80fd\u5b58\u5728\u8bef\u5dee\u3002\u5b9e\u9645\u5de5\u7a0b\u5e94\u7528\u524d\u52a1\u5fc5\u518d\u6b21\u6838\u5b9e\u6570\u636e\uff0c\u5f00\u53d1\u8005\u4e0d\u627f\u62c5\u4f7f\u7528\u8d23\u4efb\u3002"}
          </div>
          <div style={{color:C.textL,fontSize:9,marginBottom:4}}>
            {"\u672c\u5de5\u5177\u9700\u8981\u8054\u7f51\u4f7f\u7528\uff0c\u79bb\u7f51\u65f6\u65e0\u6cd5\u8ba1\u7b97 | Internet connection required for calculations"}
          </div>
          <div style={{color:C.textL,fontSize:9}}>
            {"\u6b22\u8fce\u63d0\u51fa\u5efa\u8bae\u548c\u53cd\u9988Bug | Feedback & bug reports welcome: "}
            <span style={{color:C.pri}}>chemcalc@outlook.com</span>
          </div>
        </div>
      </div>
    </div>
  );
}

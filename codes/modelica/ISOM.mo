within ;
package ISOM

  
  model hysys_comparision
     constant Chemsep_Database.Npentane comp1;
     constant Chemsep_Database.Isopentane comp2;
     constant Chemsep_Database.Nhexane comp3;
     constant Chemsep_Database.Twomethylpentane comp4;
     constant Chemsep_Database.Threemethylpentane comp5;
     constant Chemsep_Database.TwoTwodimethylbutane comp6;
     constant Chemsep_Database.TwoThreedimethylbutane comp7;
     constant Chemsep_Database.Cyclohexane comp8;
     constant Chemsep_Database.Benzene comp9;
     constant Chemsep_Database.Hydrogen comp10;
     constant Chemsep_Database.Methylcyclopentane comp11;
     constant Chemsep_Database.Cyclopentane comp12;
     constant Chemsep_Database.Nbutane comp13;
     constant Chemsep_Database.Isobutane comp14;
     
     constant Integer n = 14 "no. of components";
     constant Integer rxns = 38 "no. of reactions";
     
     constant Chemsep_Database.General_Properties comp[n]={comp1, comp2, comp3, comp4, comp5, comp6, comp7, comp8, comp9, comp10, comp11, comp12, comp13, comp14} "comp contains all the components data ";
     constant Real Kij[n,n]={{0, 0.06, 0.00118, 0, 0, 0, 0.0037, 0.0189, 0, 0, 0, 0,0.0174,0},{0.06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0056, 0},{0.00118, 0, 0, 0, 0, 0, 0, -0.003, 0.0089, -0.03, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0,0 ,0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0.037, 0,-0.003, 0, 0, 0, 0, 0, 0.0126, 0, 0, 0, 0, 0},{0.0189, 0, 0.0089, 0, 0, 0, 0, 0.0126, 0, 0, 0, 0, 0, 0},{0, 0, -0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.397, 0},{0, 0, 0, 0, 0, 0, 0,0 ,0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0,0 ,0, 0, 0, 0, 0, 0},{0.0174, -0.056, 0, 0, 0, 0, 0, 0, 0, -0.397, 0, 0, 0, -0.0004},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0004, 0}} "Binary Interaction coefficients";
     constant Real Hf0[n](each unit = "J/mol")={-146711.6263, -153649.33, -167200, -174300, -171600, -185600, -177800, -124600, 82900, 0, -106000, -77341.6263, -125600, -134200}"Heat of formation of nC5 and iC5 respectively at standard conditions";  
     constant Real Z_0[n] = {0.1491, 0.15, 0.168, 0.17, 0.167, 0.1617, 0.161, 0.1334, 0.1123, 1.01, 0.1394, 0.1146, 0.1231, 0.1269}"compressiblity factor at standard state";
     constant Real w[n] = {0.25, 0.227, 0.297, 0.278, 0.273, 0.233, 0.248, 0.211, 0.209, -0.21599, 0.227, 0.011, 0.199, 0.186} "ascentric factor";
     
     parameter Real K0[rxns](each unit="1/hr") = {4.76452E+18,7.21139E+18,33730214071.0,41497229004.0,1.04236E+21,2.93778E+21,1.34282E+19,1.34282E+19,2.50046E+14,7.04724E+15,99544918.98,8091314788.0,4.14972E+13,5.72821E+14,3.87275E+12,3.87275E+13,2177805.536,7379366.791,2.5587E+29,8.67E+25, 2.50046E+26, 8.67E+27, 1264791.964, 547040.0177, 6.42716E+09, 5861639397.0, 1.28239E+16, 1.77019E+08, 6.28086E+17, 1.31226E+10, 5.4704E+16, 8472646.706, 1.98618E+14, 867000.000, 6.73006E+20, 6.13789E+13, 4.44651E+15, 2.80556E+14}"Pre exponential factor";
     parameter Real E[rxns](each unit = "J/mol")={148.93,154.28,143.17,151.41,150.98,155.92,152.96,149.95,127.28,139.07,64.5,77.06,146.14,160.28,98.28,105.4,3.51,4.79,180.2,400.43, 187.05, 300.79, 51.08, 341.89, 129.75, 88.64,135.45, 129.29, 154.54, 98.63, 150.29, 102.35, 168.13, 90.7, 177.32, 222.9, 59.8, 54.08}"activation energy";
      
     parameter Real Fi(unit="mol/hr") = 1000*1000;
     parameter Real yi[n] = {0.04345,0.04345,0.04345,0.04345,0.04345,0.04345,0.04345,0.04345,0.04345,0.4348,0.04345, 0.04345,0.04345,0.04345} "inlet mole fraction ";
     parameter Real Ca1[n]=Fi*yi/1148;
     
     constant Integer reac1[rxns] = {1, 2, 3, 4, 3, 5, 4, 5, 4, 6, 4, 7, 5, 6, 5, 7, 6, 7, 9, 8, 9, 11, 8, 11, 3, 8, 4, 11, 5, 11, 6, 11, 7, 11, 12, 1, 13, 14};
     constant Integer prod1[rxns] = {2, 1, 4, 3, 5, 3, 5, 4, 6, 4, 7, 4, 6, 5, 7, 5, 7, 6, 8, 9, 11, 9, 11, 8, 8, 3, 11, 4, 11, 5, 11, 6, 11, 7, 1, 12, 14, 13};
     
     parameter Real P(unit="Pa")=3.2e+6 "inlet stream pressure";
     parameter Real Ti (unit="K")= 200+273.15"inlet temperature";
     
     parameter Real mi(each unit = "gm/hr") = Fi*yi*comp.MW "mass flowrate";
     parameter Real M(unit = "gm/hr") = sum(mi)"total mass flowrate";
     parameter Real xi = mi/M "initial wt fraction";
   
     constant Real T0 (unit="K") = 298.15"Standard Temp";
     constant Real Pi = 3.141592654;
     constant Real R(unit = "J/mol-K")=8.314;
     constant Real ACS(unit = "m2")=11.89/11.9;
   
     Real K[rxns](each unit="mol/hr-m3reactor");
     Real r[rxns](each unit="mol/hr-m3")"rate";
     Real Ca[n](each unit="mol/m3") "concetration";
     Real T(unit="K")"temperature";
     Real Q[rxns](each unit="J/mol")"heat evolved from reaction";
     Real denm_new;
     
     Real y_i[n]"mole fraction";
     Real x_i[n]"wt fraction";
     Real m_i[n](each unit = "g/hr")"individual mass flowrate";
     Real f_i[n](each unit="mol/hr")"molar flow rate";
     Real F(unit="mol/hr")"overall flowrate";
     Real rhom(unit = "mole/m3")"molar density";  
     Real S(unit="m3/hr")"overall Vol Flowrate";
     
     Real Cpig[n](each unit ="J/mol-K")"heat capacity of ideal components";
     Real Cp0ig[n](each unit ="J/mol-K");
     Real Cp_T1[n](each unit = "J/mol-K");
     Real Cp_T2[n](each unit = "J/mol-K");
     
     Real delH_ig[n](each unit ="J/mol");
     Real delH[rxns](each unit = "J/mol");
     Real delH_res_1[n](each unit = "J/mol");
     Real delHf0[rxns](each unit = "J/mol");
      
     Real Tr[n], Pr[n], V[n], a[n], b[n], c[n], Z[n](each start = 1), A[n], B[n], Coeff[n,4];
     Real am, bm, Zm, Am, Bm, Coeffm[4], Vm, Cpigm, Cpresm, Cpm, dadt[n], dadt_m;
     Real a_0[n], b_0[n], dadt_0[n], V_0[n];
     Real r_i[rxns], Cpm_kg, T_in_c(unit="C");
    
   
   initial equation
     Ca=Ca1;
     T=Ti;
   
   equation
   // reduced pressure and reduced temperature
     Pr = P./comp.Pc;
     Tr = T./comp.Tc;
     
   //calc of mole-fraction
     y_i = Ca/sum(Ca);
   //calc of wt-fraction
     x_i = y_i.*comp.MW/sum(y_i.*comp.MW);
   //calc of individual mass flowrate
     m_i = x_i*M;
   //calc of individual molar flow rate
     f_i = m_i./comp.MW;
   //calc of overall flowrate
     F = sum(f_i);
   //calc of molar density of stream
     rhom = P/(R*T*Zm);
   //calc of overall volumetric flowrate
     S = (1/rhom)*F;
     
   //calculating rate constant
     K = K0.*exp(-E*1000./(R*T));
   
   // calculating rate expression
     for i in 1:rxns loop
       if i<=18 then
         r[i] = K[i] * Ca[reac1[i]];
       end if;
     end for;
       r[19] = K[19] * Ca[9] * Ca[10];
       r[20] = K[20]* Ca[8];
       r[21] = K[21] * Ca[9] * Ca[10];
       r[22] = K[22] * Ca[11];
       r[23] = K[23] * Ca[8];
       r[24] = K[24] * Ca[11];
       r[25] = K[25] * Ca[3];
       r[26] = K[26] * Ca[8] * Ca[10];
       r[27] = K[27] * Ca[4];
       r[28] = K[28] * Ca[11] * Ca[10];
       r[29] = K[29] * Ca[5];
       r[30] = K[30] * Ca[11] * Ca[10];
       r[31] = K[31] * Ca[6];
       r[32] = K[32] * Ca[11] * Ca[10];
       r[33] = K[33] * Ca[7];
       r[34] = K[34] * Ca[11] * Ca[10];
       r[35] = K[35] * Ca[12] * Ca[10];
       r[36] = K[36] * Ca[1];
       r[37] = K[37] * Ca[13];
       r[38] = K[38] * Ca[14];
       
   
   // coeff calculation in the equation a*Z^3 + b*Z^2 + c*Z + d =0 and vanderwaals constant
     (Coeff, Coeffm,a, b, c, am, bm, A, B, Am, Bm) = compressiblity(P, comp.Pc, Pr, T, comp.Tc, Tr, Kij, w, y_i, n);
     ( , , a_0, b_0)  = compressiblity(P, comp.Pc, Pr, T0, comp.Tc, T0./comp.Tc, Kij, w, y_i, n);
   
   // compressiblity factor calculation
     for i in 1:n loop
       Coeff[i,1]*Z[i]^3 + Coeff[i,2]*Z[i]^2 + Coeff[i,3]*Z[i] + Coeff[i,4] = 0;
     end for;
     Coeffm[1]*Zm^3 + Coeffm[2]*Zm^2 + Coeffm[3]*Zm + Coeffm[4] = 0;
   
   // molar volume calculation
     V = (R*T/P)*Z;
     V_0 = (R*T/P)*Z_0;
     Vm = Zm*R*T/P;
     
  //   Modelica.Utilities.Streams.print(String(V_0));
     
   
   // cp values of ideal gas at different temp between T0 and T
     Cpig = Functions.VapCpId(comp.VapCp, T);
     Cp0ig = Functions.VapCpId(comp.VapCp, T0);
     Cp_T1 = Functions.VapCpId(comp.VapCp, T0+(T-T0)/3);
     Cp_T2 = Functions.VapCpId(comp.VapCp, T0+(T-T0)/6);
   // ideal cp for the entire stream
       Cpigm = sum(y_i .* Cpig);
    
   // cp residual calc
     (Cpresm, dadt_m) = Cp_res_m(P,comp.Pc,Vm,T,comp.Tc,a,b,c,am,bm,Kij,y_i,n);
   // cp real calc
       Cpm = Cpigm + Cpresm;
     
   // der(a) calc using Cp_res function which will be used in the delH_res calc
     for i in 1:n loop
         ( , , , ,dadt[i]) = Cp_res(comp[i].Pc, V[i], T, comp[i].Tc, Tr[i], w[i]);
         (  , , , ,dadt_0[i]) = Cp_res(comp[i].Pc, V_0[i], T0, comp[i].Tc, T0/comp[i].Tc, w[i]);
     end for;
   
   
     for i in 1:n loop
         //"calculating the heat required to cool reactants ideal"
         delH_ig[i] = (1/3)*(T-T0)*(Cpig[i]+Cp0ig[i]+4*(Cp_T1[i])+2*(Cp_T2[i]));
         //calc of residual H to cool reactants
         delH_res_1[i] = enthalpy_resid(a[i] ,b[i] ,Z[i] , dadt[i] ,P ,T) - enthalpy_resid(a_0[i], b_0[i], Z_0[i], dadt_0[i], P, T0);
     end for;
     
     
     for i in 1:rxns loop
       if i <=18 then
         // calc overall heat for cooling  of reactants and heating of products for each reaction
         delH[i] = delH_ig[prod1[i]]-delH_ig[reac1[i]] + delH_res_1[prod1[i]]-delH_res_1[reac1[i]];
       end if;
       //calc of heat of reaction at standard state for all the reactions
       delHf0[i] = Hf0[prod1[i]]-Hf0[reac1[i]];
     end for;
       delH[19] = (delH_ig[8] + delH_res_1[8]) - (3*(delH_ig[10] + delH_res_1[10]) + delH_ig[9] + delH_res_1[9] );
       delH[20] = -delH[19];
       delH[21] = (delH_ig[11] + delH_res_1[11]) - (3*(delH_ig[10] + delH_res_1[10]) + delH_ig[9] + delH_res_1[9] );
       delH[22] = -delH[21];
       delH[23] = delH_ig[11]-delH_ig[8] + delH_res_1[11]-delH_res_1[8];
       delH[24] = -delH[23];
       delH[25] = (delH_ig[10] + delH_res_1[10] + delH_ig[8] + delH_res_1[8]) - ( delH_ig[3] + delH_res_1[3] );
       delH[26] = -delH[25];
       delH[27] = (delH_ig[10] + delH_res_1[10] + delH_ig[11] + delH_res_1[11]) - ( delH_ig[4] + delH_res_1[4] );
       delH[28] = -delH[27];
       delH[29] = (delH_ig[10] + delH_res_1[10] + delH_ig[11] + delH_res_1[11]) - ( delH_ig[5] + delH_res_1[5] );
       delH[30] = -delH[29];
       delH[31] = (delH_ig[10] + delH_res_1[10] + delH_ig[11] + delH_res_1[11]) - ( delH_ig[6] + delH_res_1[6] );
       delH[32] = -delH[31];
       delH[33] = (delH_ig[10] + delH_res_1[10] + delH_ig[11] + delH_res_1[11]) - ( delH_ig[7] + delH_res_1[7] );
       delH[34] = -delH[33];    
       delH[35] = ( delH_ig[1] + delH_res_1[1] ) - (delH_ig[10] + delH_res_1[10] + delH_ig[12] + delH_res_1[12]);
       delH[36] = -delH[35];
       delH[37] = (delH_ig[14]+delH_res_1[14])-(delH_ig[13]+delH_res_1[13]);
       delH[38] = -delH[37];
  
      
   //heat evolved from each reaction
     Q = delHf0 + delH;
   //denominator in the enthalpy balance
     denm_new = ((F/S)*Cpm);
   
   //component balance
     der(Ca[1]) = ACS*(1/S)*(-r[1] + r[2] + r[35] - r[36]);
     der(Ca[2]) = ACS*(1/S)*(r[1] - r[2]);
     der(Ca[3]) = ACS*(1/S)*(-r[3] + r[4] - r[5] + r[6] - r[25] + r[26]);
     der(Ca[4]) = ACS*(1/S)*(r[3] - r[4] -r[7] +r[8] -r[9] + r[10] + r[11] - r[12] );
     der(Ca[5]) = ACS*(1/S)*(r[5] - r[6] + r[7] - r[8] - r[13] + r[14] - r[15] + r[16] - r[27] + r[28] - r[29] + r[30]);
     der(Ca[6]) = ACS*(1/S)*(-r[10] + r[9] + r[13] -r[14] - r[17] + r[18] -r[31] + r[32]);
     der(Ca[7]) = ACS*(1/S)*(-r[12] + r[11] + r[15] - r[16] + r[17] - r[18] -r[33] +r[34]);
     der(Ca[8]) = ACS*(1/S)*(-r[20] + r[19] - r[23] + r[24] + r[25] - r[26]);
     der(Ca[9]) = ACS*(1/S)*(r[20] - r[19] - r[21] + r[22]);
     der(Ca[10]) = ACS*(1/S)*((r[20] - r[19] - r[21] + r[22])*3  + r[25] - r[26] + r[27] - r[28] + r[29] - r[30] +r[31] - r[32] +r[33] -r[34] -r[35] + r[36]);
     der(Ca[11]) = ACS*(1/S)*(r[21] - r[22] + r[23] - r[24]  + r[27] - r[28] + r[29] - r[30] +r[31] - r[32] +r[33] -r[34]);
     der(Ca[12]) = ACS*(1/S)*(-r[35] + r[36]);
     der(Ca[13]) = ACS*(1/S)*(-r[37] + r[38]);
     der(Ca[14]) = ACS*(1/S)*(+r[37] - r[38]);
     
     
     der(T) = -ACS*(1/S)*(sum(Q.*r)/(denm_new));
     
   //unit conversions
     T_in_c = T - 273.15;
     Cpm_kg = Cpm/sum(y_i.*comp.MW);
     r_i = r*(1/3.6e+6);//converting mol/m3-hr to Kmol/m3-Sec
     
   end hysys_comparision;

function enthalpy_resid

input Real am,bm,Zm,dadT,P,T;
output Real Hr;
protected
constant Real R=8.314,
V0 = R* 298.15/101325,
uu=2,
ww=-1;
Real B,DAres,DSres, DHres;

algorithm
B := bm*P/(R*T);

DAres := am / (bm * (uu ^ 2 - 4 * ww) ^ 0.5) * log((2 * Zm + B * (uu - (uu ^ 2 - 4 * ww) ^ 0.5)) / (2 * Zm + B * (uu + (uu ^ 2 - 4 * ww) ^ 0.5))) - R * T * log((Zm - B) / Zm) - R * T * log(Zm);

DSres := R * log((Zm - B) / Zm) + R * log(Zm) - 1 / (8 ^ 0.5 * bm) * dadT * log((2 * Zm + B * (2 - 8 ^ 0.5)) / (2 * Zm + B * (2 + 8 ^ 0.5)));

Hr := DAres + T * (DSres) + R * T * (Zm - 1);

end enthalpy_resid;

function Cp_res_m
input Real P,Pc[:],Vm,T,Tc[:],a[:],b[:],c[:],am,bm,Kij[:,:],y[:],n;
output Real Cpres, dadt, d2adt, dpdt, d2pdt, dvdt;
protected
constant Real R = 8.314;
Real aux1, aux2;

algorithm

aux1 := -R/2 * (0.45724/T)^0.5;
aux2 :=0;

for i in 1:n loop
  for j in 1:n loop
    aux2 := aux2 + y[i]*y[j]*(1 - Kij[i,j])*(c[j]*(a[i]*Tc[j]/Pc[j])^0.5 + c[i]*(a[j]*Tc[i]/Pc[i])^0.5);
  end for;
end for;


dadt := aux1*aux2;
d2adt := ((R/4) * ((0.45724/T)^0.5) * (1/T))*aux2;
dpdt := R/(Vm-bm) - dadt/(Vm*(Vm+bm)+bm*(Vm-bm));
d2pdt := -d2adt/(Vm*(Vm+bm)+bm*(Vm-bm)); 
dvdt := dpdt/(R*T/(Vm-bm)^2-am*(2*Vm + 2*bm)/(Vm*(Vm+bm)+bm*(Vm-bm))^2);
Cpres := -R + T*dpdt*dvdt-T*d2adt/(8^0.5*bm)*log((Vm+(1-2^0.5)*bm)/(Vm+(1+2^0.5)*bm));
end Cp_res_m;

function Cp_res
input Real Pc,V,T,Tc,Tr,w;
output Real dpdt,d2p_dt2,dpdv,dvdt,dera,der2a,Cpres;
protected
constant Real coeff1 = (0.37464 + 1.54226*w-0.26992*w^2), coeff2=0.45724*(((R*Tc)^2)/Pc),coeff = coeff1*coeff2,
R=8.314;
Real a, b;
 
algorithm
a := ((1+coeff1*(1-Tr^0.5))^2)*coeff2;
b := 0.0778*(R*Tc)/Pc; 
dpdv := -R*T/(V-b)^2-(a*(2*V+2*b)/(V*(V+b)+b*(V-b))^2);
dera := (-coeff1*coeff2/((Tc*T)^0.5))*(1+coeff1*(1-Tr^0.5));
der2a := 0.5*((coeff1*coeff2*(1+coeff1*(1-Tr^0.5))/((Tc^0.5)*T^1.5))+(coeff1^2*coeff2/(Tc*T)));
dpdt := R/(V-b) - dera/(V*(V+b)+b*(V-b));
d2p_dt2 := -der2a/(V*(V+b)+b*(V-b)); 
dvdt := dpdt/(R*T/(V-b)^2-a*(2*V+2*b)/(V*(V+b)+b*(V-b))^2);
Cpres :=-R+T*dpdt*dvdt-T*der2a/(8^0.5*b)*log((V+(1-2^0.5)*b)/(V+(1+2^0.5)*b));
end Cp_res;

function compressiblity
input Real P, Pc[:], Pr[:], T ,Tc[:] ,Tr[:] ,Kij[:,:] ,w[:] ,y_i[:] ;
input Integer n;
output Real Coeff[n,4], Coeffm[4] ,a[n] ,b[n] ,c[n], am, bm ,A[n] , B[n] ,Am, Bm;

protected
constant Real R = 8.314;

algorithm
c := (0.37464 .+ 1.54226*w .- 0.26992*w.^2);
a := (1 .+ c .* (1 .- Tr.^0.5)).^2 .* (0.45724 .* (((R*Tc).^2) ./ Pc));
b := 0.0778*(R*Tc)./Pc;

A := a * (P/(R*T)^2);
B := b * (P/(R*T));


Coeff[:,1] := ones(n);
Coeff[:,2] := B .- 1;
Coeff[:,3] := A .- 3*(B.^2) .- 2*B ;
Coeff[:,4] := B.^3 + B.^2 .- A.*B;

am := 0.0;
bm := 0.0;
for i in 1:n loop
    for j in 1:n loop
      am := am + y_i[i]*y_i[j]*sqrt(a[i]*a[j])*(1-Kij[i,j]);
    end for;
   bm := bm + y_i[i]*b[i];
end for;

Am := am * (P/(R*T)^2);
Bm := bm * (P/(R*T));

Coeffm[1] := 1;
Coeffm[2] := Bm-1;
Coeffm[3] := Am - 3*Bm^2 - 2*Bm;
Coeffm[4] := Bm^3 + Bm^2 - Am*Bm;

end compressiblity;



end ISOM;

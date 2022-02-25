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
        
        constant Integer n = 10 "no. of components";
        constant Integer rxns = 16 "no. of reactions";
        
        constant Chemsep_Database.General_Properties comp[n]={comp1, comp2, comp3, comp4, comp5, comp6, comp7, comp8, comp9, comp10} "comp contains all the components data ";
        constant Real Kij[n,n]={{0, 0.06, 0.00118, 0, 0, 0, 0.0037, 0.0189, 0, 0},{0.06, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0.00118, 0, 0, 0, 0, 0, 0, -0.003, 0.0089, -0.03},{0, 0, 0, 0, 0, 0, 0,0 ,0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0,0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0.037, 0,-0.003, 0, 0, 0, 0, 0, 0.0126, 0},{0.0189, 0, 0.0089, 0, 0, 0, 0, 0.0126, 0, 0},{0, 0, -0.03, 0, 0, 0, 0, 0, 0, 0}} "Binary Interaction coefficients";
        constant Real Hf0[n](each unit = "J/mol")={-146711.6263, -153649.33, -167200, -174300, -171600, -185600, -177800, -124600, 82900, 0}"Heat of formation of nC5 and iC5 respectively at standard conditions";   
        parameter Real K0[rxns](each unit="1/hr") = { 4.76452E+18, 7.21139E+18, 33730214071.0, 41497229004.0, 1.04236E+21, 2.93778E+21, 1.34282E+19, 1.34282E+19, 2.50046E+14, 7.04724E+15, 99544918.98, 8091314788.0, 4.14972E+13, 5.72821E+14 ,6.42716E+9, 5861639397.0}"Pre exponential factor";
        parameter Real E[rxns](each unit = "J/mol")={148930, 154280, 143170, 151410.0, 150980, 155920, 152960, 149950, 127280, 139070, 64500, 77060, 146140, 160280 ,129750, 88640}"activation energy";
        constant Real Z_0[n] = {0.1491, 0.15, 0.168, 0.17, 0.167, 0.1617, 0.161, 0.1334, 0.1123, 1.01}"compressiblity factor at standard state";
        constant Real w[n] = {0.25, 0.227, 0.297, 0.278, 0.273, 0.233, 0.248, 0.211, 0.209, -0.21599} "ascentric factor";
        
        parameter Real P(unit="Pa")=3.2e+6 "inlet stream pressure";
        parameter Real Ti (unit="K")= 200+273.15"inlet temperature";
         
        parameter Real Fi(unit="mol/hr") = 1000*1000;
        parameter Real yi[n] = {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.55} "inlet mole fraction ";
        parameter Real Ca1[n]=Fi*yi/1148;
        
        constant Integer reac1[rxns] = {1, 2, 3, 4, 3, 5, 4, 5, 4, 6, 4, 7, 5, 6, 9, 8};
        constant Integer prod1[rxns] = {2, 1, 4, 3, 5, 3, 5, 4, 6, 4, 7, 4, 6, 5, 8, 9};
        
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
        K = K0.*exp(-E./(R*T));
      
      // calculating rate expression
        for i in 1:rxns loop
          if i<=14 then
            r[i] = K[i] * Ca[reac1[i]];
          end if;
        end for;
          r[15] = K[15] * Ca[9] * Ca[10];
          r[16] = K[16]* Ca[8];
      
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
          if i <=14 then
            // calc overall heat for cooling  of reactants and heating of products for each reaction
            delH[i] = delH_ig[prod1[i]]-delH_ig[reac1[i]] + delH_res_1[prod1[i]]-delH_res_1[reac1[i]];
          end if;
          //calc of heat of reaction at standard state for all the reactions
          delHf0[i] = Hf0[prod1[i]]-Hf0[reac1[i]];
        end for;
          delH[15] = (delH_ig[8] + delH_res_1[8]) - (3*(delH_ig[10] + delH_res_1[10]) + delH_ig[9] + delH_res_1[9] );
          delH[16] = -delH[15];
         
      //heat evolved from each reaction
        Q = delHf0 + delH;
      //denominator in the enthalpy balance
        denm_new = ((F/S)*Cpm);
      
      //component balance
        der(Ca[1]) = ACS*(1/S)*(-r[1] + r[2]);
        der(Ca[2]) =  ACS*(1/S)*(r[1] - r[2]);
        der(Ca[3]) =  ACS*(1/S)*(-r[3] + r[4] - r[5] + r[6]);
        der(Ca[4]) = ACS*(1/S)*(r[3] - r[4] -r[7] +r[8] -r[9] + r[10] +r[11] - r[12] );
        der(Ca[5]) = ACS*(1/S)*(r[5] - r[6] + r[7] - r[8] -r[13] + r[14]);
        der(Ca[6]) = ACS*(1/S)*(-r[10] + r[9] + r[13] -r[14]);
        der(Ca[7]) = ACS*(1/S)*(-r[12] + r[11]);
        der(Ca[8]) = ACS*(1/S)*(-r[16] + r[15]);
        der(Ca[9]) = ACS*(1/S)*(r[16] - r[15]);
        der(Ca[10]) = ACS*(1/S)*(r[16] - r[15])*3;
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

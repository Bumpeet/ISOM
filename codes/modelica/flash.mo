package flash
  model imp
  parameter Chemsep_Database.Benzene comp1;
  parameter Chemsep_Database.Toluene comp2;
    parameter Integer n = 2;
    parameter Chemsep_Database.General_Properties comp[n] = {comp1, comp2};
    parameter Real z[n] = {0.5, 0.5};
    parameter Real T = 100 + 273.15;
    parameter Real P = 1.1e+5;
    parameter Real F(unit = "mol/hr") = 100;
    parameter Real w[n] = comp.AF;
    Real Pc[n], Tc[n];
    Real a[n], b[n];
    Real Coeffm[2, 4], Zm[2, 3], am[2], bm[2], Am[2], Bm[2], aij[n, n];
    Real phi[n, 2];
    Real y[n](each min=0, each max=1), x[n](each min=0, each max=1);
    Real L(start = 60, min=0, max=100), V(start = 40, min=0, max=100);
    Real Ax[n],Ay[n];
    Real k[2], alpha[2];
  equation
    Pc = comp.Pc;
    Tc = comp.Tc;
    F = L + V;
    
    for i in 1:n loop
      F * z[i] = L * x[i] + V * y[i];
      y[i] / x[i] = phi[i, 1] / phi[i, 2];
      Ax[i] = sum(aij[i,:].*x);
      Ay[i] = sum(aij[i,:].*y);
      log(phi[i, 1]) = b[i] / bm[1] * (Zm[1, 1] - 1) - log(Zm[1, 1] - Bm[1]) - (Am[1] / (2 * sqrt(2)* Bm[1])) * (2 * Ax[i] / Am[1] - b[i] / bm[1]) * log((Zm[1, 1] + 2.414 * Bm[1]) / (Zm[1, 1] - 0.414 * Bm[1]));
      log(phi[i, 2]) = b[i] / bm[2] * (Zm[2, 3] - 1) - log(Zm[2, 3] - Bm[2]) - (Am[2] / (2 * sqrt(2)* Bm[2])) * (2 * Ay[i]/ Am[2] - b[i] / bm[2]) * log((Zm[2, 3] + 2.414 * Bm[2]) / (Zm[2, 3] - 0.414 * Bm[2]));
    end for;
    sum(y) = sum(x);
    (a, b, k, alpha) = coeff(T, Tc, P, Pc, w);
    (Coeffm[1, :], am[1], bm[1], Am[1], Bm[1], aij) = coeffm(T, P, x, a, b);//coeff for liquid phase
    (Coeffm[2, :], am[2], bm[2], Am[2], Bm[2]) = coeffm(T, P, y, a, b);//coeff for gaseous phase
    Zm = sort(2,cubic_root(Coeffm));//we'll get 3 roots for each phase
    
  end imp;

  function coeff
    input Real T, Tc[:], P, Pc[:], w[:];
    output Real a[2], b[2], k[2], alpha[2];
  protected
  constant Integer n = 2;
    
    Real A[n], B[n];
    Real Tr[n];
    constant Real R = 8.314;
  algorithm
    Tr := T./Tc;
    for i in 1:n loop
    k[i] := 0.37464 + 1.54226 * w[i]- 0.26992 * (w[i]^ 2);
    alpha[i] := (1 + k[i] * ((1 - Tr[i] ^ 0.5)) )^ 2;
    a[i] :=  alpha[i] * 0.45724 * ((R * Tc[i]) ^ 2) / Pc[i];
    b[i] := 0.0778 * (R * Tc[i]) / Pc[i];
  
    end for;
  //  a[1]:=2.6;
  //  a[2]:=3.59;
  end coeff;

  function cubic_root
    input Real coeff[2,4];
    output Real soln[2,3];
  protected
    constant Integer n = size(coeff,1);
    Real a, b, c;
    Real Q, R, theta;
    Real S, T, M;
  algorithm
    for i in 1:n loop
      a := coeff[i,2] / coeff[i,1];
      b := coeff[i,3] / coeff[i,1];
      c := coeff[i,4] / coeff[i,1];
      Q := (a ^ 2 - 3 * b) / 9;
      R := (2 * a ^ 3 - 9 * a * b + 27 * c) / 54;
      M := R ^ 2 - Q ^ 3;
      if M <= 0 then
        theta := Modelica.Math.acos(R / ((Q ^ 3) ^ 0.5));
        soln[i,1] := (-2 * Q ^ 0.5 * Modelica.Math.cos(theta / 3)) - a / 3;
        soln[i,2] := (-2 * Q ^ 0.5 * Modelica.Math.cos((theta + 2 * 3.14) / 3)) - a / 3;
        soln[i,3] := (-2 * Q ^ 0.5 * Modelica.Math.cos((theta - 2 * 3.14) / 3)) - a / 3;
      else
        S := ((-R) + M ^ 0.5) ^ (1 / 3);
        T := ((-R) - M ^ 0.5) ^ (1 / 3);
        soln[i,1] := S + T - a / 3;
        soln[i,2] := S + T - a / 3;
        soln[i,3] := S + T - a / 3;
      end if;
     end for;
  end cubic_root;

  function coeffm
    input Real T, P, y_i[:];
    input Real a[:], b[:];
    output Real Coeffm[4];
    output Real am, bm, Am, Bm, aij[2,2];
  protected
    constant Real R = 8.314;
    constant Integer n = 2;
  algorithm
    am := 0.0;
    bm := 0.0;
    for i in 1:n loop
      for j in 1:n loop
        aij[i, j] := sqrt(a[i] * a[j]);
        am := am + y_i[i] * y_i[j] * aij[i,j];
      end for;
      bm := bm + y_i[i] * b[i];
    end for;
    Am := am * (P / (R * T) ^ 2);
    Bm := bm * (P / (R * T));
    aij := aij*(P / (R * T) ^ 2);
    Coeffm[1] := 1.0; 
    Coeffm[2] := Bm - 1;
    Coeffm[3] := Am - 3 * Bm ^ 2 - 2 * Bm;
    Coeffm[4] := Bm ^ 3 + Bm ^ 2 - Am * Bm;
  end coeffm;

  function sort
    input Integer n;
    input Real x[:, 3];
    output Real y[2, 3];
  algorithm
    for i in 1:n loop
      if x[i, 1] <= x[i, 2] then
        if x[i, 2] <= x[i, 3] then
          y[i, :] := x[i, :];
        else
          y[i, 1] := x[i, 1];
          y[i, 2] := x[i, 3];
          y[i, 3] := x[i, 2];
        end if;
      elseif x[i, 1] > x[i, 2] then
        if x[i, 2] >= x[i, 3] then
          y[i, 1] := x[i, 3];
          y[i, 2] := x[i, 2];
          y[i, 3] := x[i, 1];
        else
          y[i, 1] := x[i, 2];
        end if;
        if x[i, 1] < x[i, 3] then
          y[i, 2] := x[i, 1];
          y[i, 3] := x[i, 3];
        else
          y[i, 2] := x[i, 3];
          y[i, 3] := x[i, 1];
        end if;
      end if;
    end for;
  end sort;

  model flash_algo
    parameter Chemsep_Database.Benzene comp1;
    parameter Chemsep_Database.Toluene comp2;
    parameter Integer n = 2;
    parameter Chemsep_Database.General_Properties comp[n] = {comp1, comp2};
    parameter Real z[n] = {0.4, 0.6};
    parameter Real T = 273 + 100;
    parameter Real P = 101325;
    parameter Real F(unit = "mol/hr") = 100;
    parameter Real w[n] = comp.AF;
    
    constant Real eps = 1e-6;
    
    Real Pc[n], Tc[n];
    Real Coeff[n, 4], a[n], b[n];
    Real CoeffmL[4], CoeffmV[4], Zm[2,3], am[2], bm[2], Am[2], Bm[2], aij[n, n];
    Real phi[n, 2];
    Real y[n](start = {0.25, 0.75}), x[n](start = {0.75, 0.25});
    Real psi(start=0.4), K[2], fun_psi;
  equation
  Pc = comp.Pc;
  Tc = comp.Tc;
  
  algorithm
    (Coeff, a, b) := coeff(T, Tc, P, Pc, w);
    (CoeffmL, am[1], bm[1], Am[1], Bm[1], aij) := coeffm(T, P, x, a, b);//coeff for liquid phase
    (CoeffmV, am[2], bm[2], Am[2], Bm[2]) := coeffm(T, P, y, a, b);//coeff for gaseous phase
  //  ZmL := cubic_root(CoeffmL);//we'll get 3 roots for each phase
  //  ZmV := cubic_root(CoeffmV);
    Zm := cubic_root({CoeffmL, CoeffmV});
    
    for i in 1:n loop
      (phi[i, 1]) := exp(b[i] / bm[1] * (Zm[1,1] - 1) - log(Zm[1,1] - Bm[1]) - (Am[1] / (2 * 2 ^ 0.5 * Bm[1])) * (2 * sum(aij[i, :] .* x) / am[1] - b[i] / bm[1]) * log((Zm[1,1] + 2.414 * Bm[1]) / (Zm[1,1] - 0.414 * Bm[1])));
      (phi[i, 2]) := exp(b[i] / bm[2] * (Zm[2,1] - 1) - log(Zm[2,1] - Bm[2]) - (Am[2] / (2 * 2 ^ 0.5 * Bm[2])) * (2 * sum(aij[i, :] .* y) / am[2] - b[i] / bm[2]) * log((Zm[2,1] + 2.414 * Bm[2]) / (Zm[2,1] - 0.414 * Bm[2])));
      K[i] := phi[i,1]/phi[i,2];
      x[i] := z[i]/(1+(K[i]-1)*psi);
      y[i] := K[i]*x[i];
    end for;
    
    fun_psi := abs(sum(y) - sum(x));
    if fun_psi > eps then
      psi := psi*eps/(fun_psi);
    end if;
    
  end flash_algo;
end flash;

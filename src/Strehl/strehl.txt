

"""
From Mathematica Notebook

FFTShift[list_] := Block[{n, m},
   {n, m} = Quotient[Dimensions[list], 2];
   RotateRight[#, m] & /@ (Transpose[RotateRight[#, n] & /@ Transpose[list]])];

GenerateOTF[SIZE_(*px*), LAM_(*\[Mu]m*), DLAM_(*\[Mu]m*), 
   PSCALE_(*mas/px*), M1_(*m*), M2_(*m*)] := 
  Block[{Sinc, H1, H2, G, TelOTF, lam, dlam, pscale, lambda, fmax, fc,
     x, y, r, f},
   Sinc[x_] := If[Abs[x] < 10^-4, 1, Sin[x]/x]; 
   H1[f_, u_, v_] := Block[{e},
     If[Abs[v - 1] < 10^-12, e = 1, e = -1];
     (v^2/\[Pi]) ArcCos[(f/v) (1 + e (1 - u^2)/(4 f^2))]
     ];
   H2[f_, u_] := Block[{a, b},
     a = 2 f/(1 + u);
     b = (1 - u)/(2 f);
     -1 (f/\[Pi]) (1 + u) Sqrt[(1 - a^2) (1 - b^2)]
     ];
   G[f_, u_] := Which[
     f <= (1 - u)/2, u^2,
     f >= (1 + u)/2, 0,
     True, H1[f, u, 1] + H1[f, u, u] + H2[f, u]
     ];
   TelOTF[f_, u_] := (G[f, 1] + u^2 G[f/u, 1] - 2 G[f, u])/(1 - u^2);
   lam = LAM 10^-6(*m*);
   dlam = DLAM 10^-6(*m*);
   pscale = PSCALE (2 \[Pi]/360/60^2/10^3)(*rad*);
   fmax = M1 pscale SIZE/lam;
   Total[Table[
     lambda = lam - dlam (k - 5)/8;
     fc = fmax lam/lambda;
     Array[Function[{j, i},
       y = j - SIZE/2;
       x = i - SIZE/2;
       r = Sqrt[x^2 + y^2];
       f = r/fc;
       If[f < 1, 
         If[r < 0.1, 1, 
          TelOTF[f, M2/M1] Sinc[\[Pi] x/SIZE] Sinc[\[Pi] y/SIZE]], 0]/
        9],
      {SIZE, SIZE}, {0.5, 0.5}],
     {k, 1, 9}]]];
(* ESO Eclipse 5.0.0 *)


OTF2PSF[otf_] := Block[{psf},
   psf = Abs[Fourier[otf, FourierParameters -> {1, -1}]];
   psf = FFTShift[psf];
   psf];

ModelPSF[SIZE_, LAM_, DLAM_, PSCALE_, M1_, M2_] := 
  Image@OTF2PSF[GenerateOTF[SIZE, LAM, DLAM, PSCALE, M1, M2]];


"

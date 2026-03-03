(* ================================================================
   HydrogenPlot.m
   2D Probability Density Plots for Hydrogen-like Atoms

   Atomic units throughout: a0 = hbar = me = e = 1, Z = 1.
   Wavefunction: psi_{n,l,m}(r,theta,phi) = R_{n,l}(r) * Y_l^m(theta,phi)

   Quantum numbers:
     n  – principal quantum number  : 1, 2, 3, ...
     l  – orbital angular momentum  : 0, 1, ..., n-1
     m  – magnetic                  : -l, ..., l
     s  – spin projection (ms)      : +1/2 or -1/2

   Reference: Griffiths, "Introduction to Quantum Mechanics", Eq. 4.73
   ================================================================ *)


(* ----------------------------------------------------------------
   Section 1: Wavefunction Components
   ---------------------------------------------------------------- *)

(* Radial wavefunction R_{n,l}(r) in atomic units (Z=1).
   Uses the generalized Laguerre polynomial LaguerreL[k, alpha, x]. *)
RadialWF[n_Integer, l_Integer, r_?NumericQ] :=
  With[{rho = 2.0 r / n},
    Sqrt[(2.0/n)^3 * (n - l - 1)! / (2.0 n * (n + l)!)] *
      Exp[-r/n] * rho^l * LaguerreL[n - l - 1, 2 l + 1, rho]
  ]

(* Full wavefunction psi_{n,l,m}(r, theta, phi).
   SphericalHarmonicY[l, m, theta, phi] is Mathematica's built-in Y_l^m. *)
HydrogenWF[n_Integer, l_Integer, m_Integer,
           r_?NumericQ, theta_?NumericQ, phi_?NumericQ] :=
  RadialWF[n, l, r] * SphericalHarmonicY[l, m, theta, phi]

(* Probability density |psi_{n,l,m}|^2 at Cartesian coordinates (x, z).
   The cross-section is taken in the xz-plane (phi = 0).
   Note: |psi|^2 is phi-independent because |e^{i m phi}| = 1,
   so the xz cross-section is representative of any meridional plane. *)
ProbDensity[n_Integer, l_Integer, m_Integer, x_?NumericQ, z_?NumericQ] :=
  Module[{r, theta},
    r = Sqrt[x^2 + z^2];
    If[r < 1*^-10,
      0.0,
      theta = ArcCos[Clip[z/r, {-1., 1.}]];
      Abs[HydrogenWF[n, l, m, r, theta, 0.0]]^2
    ]
  ]


(* ----------------------------------------------------------------
   Section 2: HydrogenPlot — Main Plotting Function
   ----------------------------------------------------------------

   HydrogenPlot[n, l, m, s]

   Generates a DensityPlot of |psi_{n,l,m}|^2 in the xz-plane.

   Arguments:
     n : principal quantum number       (integer >= 1)
     l : orbital angular momentum       (integer, 0 <= l <= n-1)
     m : magnetic quantum number        (integer, -l <= m <= l)
     s : spin projection ms             (+1/2 or -1/2)
         The spin quantum number does not affect the spatial
         probability density; it is included in the plot label only.

   Returns: a DensityPlot graphic, or $Failed on invalid input.
   ---------------------------------------------------------------- *)

HydrogenPlot::usage =
  "HydrogenPlot[n, l, m, s] plots the probability density |psi_{n,l,m}|^2 " <>
  "in the xz-plane (atomic units, Z = 1). " <>
  "Spin s = +1/2 or -1/2 is noted in the title but does not affect |psi|^2. " <>
  "Option \"Gamma\" (default 0.3) applies a power-law compression density^gamma " <>
  "to prevent the nuclear peak from washing out orbital structure.";

(* Gamma < 1 compresses the dynamic range of the colour scale.
   0.3 works well across all n; raise toward 1.0 for a more linear scale. *)
Options[HydrogenPlot] = {"Gamma" -> 0.3};

HydrogenPlot::badQN =
  "Invalid quantum numbers: n = `1`, l = `2`, m = `3`. " <>
  "Requirements: n >= 1, 0 <= l <= n-1, |m| <= l.";

(* Orbital sublevel labels *)
$OrbitalLabels = <|0 -> "s", 1 -> "p", 2 -> "d", 3 -> "f", 4 -> "g", 5 -> "h"|>;

HydrogenPlot[n_Integer, l_Integer, m_Integer, s_, opts : OptionsPattern[]] :=
  Module[{rMax, lLabel, spinStr, titleStr, gamma},
    gamma = OptionValue["Gamma"];

    (* Validate quantum numbers *)
    If[n < 1 || l < 0 || l > n - 1 || Abs[m] > l,
      Message[HydrogenPlot::badQN, n, l, m];
      Return[$Failed]
    ];

    lLabel  = Lookup[$OrbitalLabels, l, "?"];
    spinStr = If[s === 1/2, "+\[HBar]/2", "-\[HBar]/2"];

    (* Adaptive range: covers > 99% of the radial probability for shell n *)
    rMax = 5.0 n^2;

    titleStr =
      "n = " <> ToString[n] <>
      ",  l = " <> ToString[l] <> " (" <> lLabel <> ")" <>
      ",  m = " <> ToString[m] <>
      ",  \!\(\*SubscriptBox[\(m\), \(s\)]\) = " <> spinStr;

    DensityPlot[
      ProbDensity[n, l, m, x, z]^gamma,
      {x, -rMax, rMax},
      {z, -rMax, rMax},
      PlotPoints          -> 120,
      ColorFunction       -> "SunsetColors",
      ColorFunctionScaling -> True,
      PlotLabel           -> Style[titleStr, 13, Bold],
      FrameLabel          -> {
        {Style["z / \!\(\*SubscriptBox[\(a\), \(0\)]\)", 11], None},
        {Style["x / \!\(\*SubscriptBox[\(a\), \(0\)]\)", 11], None}
      },
      Frame               -> True,
      AspectRatio         -> 1,
      ImageSize           -> 360
    ]
  ]


(* ----------------------------------------------------------------
   Section 3: QuantumNumberTable — Valid States up to n = nMax
   ----------------------------------------------------------------

   QuantumNumberTable[nMax]

   Returns a formatted Grid showing all valid combinations of
   quantum numbers (n, l, m, ms) for n = 1 to nMax, together
   with the corresponding Dirac ket |n, l, m, ms>.

   Rows are grouped and colour-coded by principal quantum number n.
   ---------------------------------------------------------------- *)

QuantumNumberTable[nMax_Integer] :=
  Module[{states = {}, rows, header, nPalette, rowColors},

    (* Collect all valid (n, l, m, ms) in standard ordering *)
    Do[
      AppendTo[states, {n, l, m, s}],
      {n, 1, nMax}, {l, 0, n - 1}, {m, -l, l}, {s, {1/2, -1/2}}
    ];

    (* Format each state as a table row *)
    rows = Map[
      Function[st,
        Module[{n, l, m, s, ll, msStr, ketStr},
          {n, l, m, s} = st;
          ll     = Lookup[$OrbitalLabels, l, "?"];
          msStr  = If[s === 1/2, "+1/2", "\[Minus]1/2"];
          ketStr = Row[{
            "\[LeftBracketingBar]",
            n, ", ", l, ", ", m, ", ",
            If[s === 1/2,
              "+\!\(\*FractionBox[\(1\), \(2\)]\)",
              "\[Minus]\!\(\*FractionBox[\(1\), \(2\)]\)"
            ],
            "\[RightAngleBracket]"
          }];
          {n, Row[{l, " (", ll, ")"}], m, msStr, ketStr}
        ]
      ],
      states
    ];

    (* Header row *)
    header = Style[#, Bold, 11, White] & /@ {
      "n", "l", "m",
      Subscript["m", "s"],
      "State \[LeftBracketingBar]n, l, m, \!\(\*SubscriptBox[\(m\), \(s\)]\)\[RightAngleBracket]"
    };

    (* Row background colours, grouped by n *)
    nPalette = <|
      1 -> Lighter[Blue,   0.72],
      2 -> Lighter[Teal,   0.72],
      3 -> Lighter[Green,  0.72],
      4 -> Lighter[Orange, 0.72]
    |>;
    rowColors = Prepend[
      Lookup[nPalette, #[[1]], White] & /@ states,
      GrayLevel[0.25]   (* header background *)
    ];

    Grid[
      Prepend[rows, header],
      Frame      -> All,
      FrameStyle -> Directive[GrayLevel[0.55], Thin],
      Alignment  -> {Center, Center},
      Spacings   -> {2, 0.65},
      Background -> {None, rowColors},
      ItemStyle  -> Directive[FontSize -> 10, FontFamily -> "Courier"]
    ]
  ]


(* ================================================================
   Section 4: Output
   ================================================================ *)

(* --- Quantum number table for n = 1 to 4 --- *)
Print[Style[
  "Hydrogen Atom \[LongDash] All Quantum Number Combinations (n \[LessEqual] 4)",
  16, Bold, Blue
]]
Print[QuantumNumberTable[4]]
Print[""]
Print[Style[
  "Total states: " <> ToString[2 * 4^2] <> "  (2n\[Squared] = " <>
  ToString[2*1^2] <> " + " <> ToString[2*2^2] <> " + " <>
  ToString[2*3^2] <> " + " <> ToString[2*4^2] <> " = 60)",
  12, Italic, GrayLevel[0.3]
]]

(* --- Example probability density plots --- *)
Print[""]
Print[Style["Probability Densities \[LeftBracketingBar]\[Psi]\[RightBracketingBar]\[Squared] \[LongDash] Selected Orbitals", 14, Bold]]

GraphicsGrid[{
  {HydrogenPlot[1, 0,  0,  1/2],   HydrogenPlot[2, 0,  0,  1/2]},
  {HydrogenPlot[2, 1,  0,  1/2],   HydrogenPlot[2, 1,  1,  1/2]},
  {HydrogenPlot[3, 2,  0,  1/2],   HydrogenPlot[3, 2,  1,  1/2]},
  {HydrogenPlot[4, 3,  0,  1/2],   HydrogenPlot[4, 3,  2,  1/2]}
},
  ImageSize -> 800,
  Frame     -> All,
  Spacings  -> {5, 5}
]

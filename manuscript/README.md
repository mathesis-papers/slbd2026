# KD-real SLD counterexample and mixed-state exclusions

Revision 2. Posing locator: Spriet–Langrenez–Brummelhuis–De Bièvre,
arXiv:2511.01996v2, §6.3 immediately after (6.23), p. 43; the all-mixed-state
question is on p. 45. Journal locator: DOI 10.3390/sym18061008, §6.3 after (168).

## Results

Theorem 1 gives faithful KD-positive dimension-nine states with KD-real SLD
and QFI `2 lambda^2 / ((2+lambda)(3+lambda))`, for real `lambda > -1`.
At `lambda=1` the QFI is `1/6`.
Corollary 2 extends the construction to every square dimension `n^2`, `n>=3`:
`4 (n-gcd(2,n)) lambda^2 / ((n-1)^2 (2+lambda) (n+lambda))`.
Corollary 3 adds rank-six and rank-three endpoints with QFI 1 and 2.
Corollaries 4–7 describe an optimal KD-real measurement and exclusion for
basis-diagonal spans, real transition matrices up to basis phases, cyclic-prime
Fourier bases and Haar-almost every pair of bases in each finite dimension.
In particular, the least counterexample dimension is between three and nine.
Proposition 8 proves the star-algebra/Jordan-map mechanism.

`coverage.json` and the PDF table map the posing-source statements to proved,
refuted-and-corrected, prior proof cited, or not covered. `metadata.json` lists
the single associated problem identifier, `slbd2026`.

## Build and exact checks

Use Tectonic **0.17.0**:

```sh
tectonic note.tex
python3 verify.py > verification.actual.json
cmp verification.actual.json verification.json
```

Alternatively, run `pdflatex note.tex` twice. References are included directly
in the TeX source, so BibTeX is unnecessary. The checker uses Python 3's standard
library and exact arithmetic over `Q(sqrt(3), i)`; it checks dimension-nine
identities, six rational parameters, the two endpoints, square orders 2, 4 and 6,
and a basis-diagonal trace obstruction. The continuous and arbitrary-dimension
statements have written proofs in the note.

## Formal coverage

The PDF links the three corresponding `.lean` files at revision
`71b4592f31a0ce564acedcbee195fbd50f30bd8c`. In
`KirkwoodDiracSLDRefutation.lean`,
`D5.S3.Quantum.Information.KirkwoodDiracSLDRefutation.SprietKirkwoodDiracSLD.result`
proves the negation of the mixed-state claim using the concrete dimension-nine
`lambda=1` witness. `FisherDiscord/Defs.lean` defines `IsSLD` and `fisherSLD`;
the latter is one quarter of the QFI used in the note.
`SLDSpectralEnergy.lean` defines the density predicate. These declarations do
not assert the parameter, endpoint, square, exclusion or algebra extensions.

## Creators and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.
Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.
Creative Commons Attribution 4.0 International, **CC BY 4.0** (`cc-by-4.0`).

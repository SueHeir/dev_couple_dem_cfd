# References

The measured bed force uses the independent MacDonald packed-bed closure through
the live DEM-CFD seam. The validation reference is Wen & Yu's minimum
fluidization correlation:

```text
Re_mf = sqrt(33.7^2 + 0.0408 Ar) - 33.7
Ar    = rho_f (rho_p - rho_f) g d^3 / mu^2
U_mf  = Re_mf mu / (rho_f d)
```

Sources:

- C.Y. Wen and Y.H. Yu, "A generalized method for predicting the minimum
  fluidization velocity", AIChE Journal 12(3):610-612, 1966.
- I.F. MacDonald, M.S. El-Sayed, K. Mow, F.A.L. Dullien, "Flow through Porous
  Media - the Ergun Equation Revisited", Industrial & Engineering Chemistry
  Fundamentals 18(3):199-208, 1979.
- S. Ergun, "Fluid flow through packed columns", Chemical Engineering Progress
  48(2):89-94, 1952.

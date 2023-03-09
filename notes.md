Notes
-----

`wtot` = is the sum of all the off-diagonal terms, I think it is written such that they add the off-diagonal terms to the potential by subtracting a diagonal bit and then adding a diagonal + off-diagonal bit (so that the net effect is just adding the off-diagonal part)  
`wxdsic` = the off-diagonal contribution to the potential  
`wrefsic` = the self-Hartree potential, plus other contributions that only appear in the K and KPZ functionals

Why are they members of nksiclib and not simply local variables?

Also something strange is going on with intent: e.g. in  `nksic_getOmattot_general` is `intent(in)` but internally it calls `nksic_potential` which generates `wtot` with `intent(out)`?!

Why do they need to be kept separate from vsic? Ultimately they're added to vsic

Do I need to worry about `ngm` vs `ngmx`?




# F* audit PR #60 @ ee3564ff

Source: src/sat/smt/arith_theory_checker.h at head (sha256 229f6a7a622535b75d3cd72692590e859ad5be1642fcb9b82ed4b112f6eecb58, expand_definitions lines 458-~480); base 167e71a6.

## Obligations (all proved)
1. O1: for k nonzero, u=k*d+m, 0<=m<|k|, the expansion c*s*u - c*s*k*d equals c*rem(u,k) with rem = m (k>0) or -m (k<0), any coefficient c. Matches the guard is_rem && int numeral && nonzero, s=-1 iff n<0.
2. O2: rem range: k>0 => 0<=rem(k; k<0 =) k(rem<=0.
3. O3: concrete witness u=5,k=-3 shows rem=-2 != mod=2, confirming the wrong-sign identity must be rejected (as PR tests expect).

## Correspondence
F* ints model integer numerals; div/mod defined by the SMT-LIB identities (Z3 rem = ite(k)=0, mod, -mod), as stated by the PR). Assumed: that Z3 rem semantics, the linearize/add row helpers are linear and correct, hash-consing of idiv terms so they cancel.

## Excluded / manual
Tests (src/test/proof_replay.cpp), Python test, docs, linearize/add plumbing, memory safety, whole-solver soundness; rem bounds (rejected, unchanged). Zero and term divisors: guard excludes them by code reading only (not mechanized).

## Reproduce
fstar.exe --smt z3-4.13.3/bin/z3 --z3version 4.13.3 --z3rlimit 10 --cache_dir cache --odir out RemExpansionAudit.fst (F* 2026.09.27, Z3 4.13.3). No admits/assumes.
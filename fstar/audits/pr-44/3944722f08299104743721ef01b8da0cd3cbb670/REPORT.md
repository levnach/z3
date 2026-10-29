# F* audit of PR #44 @ 3944722f08299104743721ef01b8da0cd3cbb670 (base f87e9c80)

Tool: F* 2026.09.27, Z3 4.13.3 (pinned independent backend). Module Z3ArithAxiomAudit.fst; no admits/assumes; no base modules imported (only F* standard library).

## Obligations (all proved)
1. o1_expand_mod: c*(mod u n) = c*u - (c*n)*(div u n) given SMT-LIB div/mod definition (n<>0, 0<=m<|n|). Matches expand_definitions mod branch.
2. o1_expand_div: k*(c*q)=c*u given k*q=u, k<>0; matches the (u / k) -> (1/k)*u branch.
3. o2_infeasible: if check_ineq holds after substituting mod = |n|-1 (c<0) or 0 (c>=0), then c*m+const <=0 (or <0 if strict) is unsatisfiable for m in [0,|n|-1]. Matches is_infeasible.
4. o3_cancel: positive multiple of s*u<0 added to the negated literal cancels all variables with nonnegative constant, so the two cannot hold together. Matches implies_nonnegative.
5. o3_abs: s in {1,-1}, s*u>=0 implies abs u = s*u. Matches the case split in is_abs_case.

## Correspondence and assumptions
F* models mathematical integers and the row arithmetic abstractly; the C++ side (linearize, obj_map handling, rational arithmetic, is_numeral, strictness flag for integer +1 adjustment in add_literal, m_pinned lifetime) was manually reviewed, not mechanized. Trivial literals true/false and equality literal handling (is_valid_literal) were read manually: sound because positive eq needs zero row, negated eq needs nonzero constant row after expansion. Dropping a true argument in clause() was manually reviewed as the q=0 literal being false.

## Exclusions
proof_cmds.cpp, tests (proof_replay.cpp, test_proof_certificate.py), docs, rule registration/plumbing, check_axiom loop, real-valued abs/mod edge cases, memory safety, and the rest of the solver. Findings: none (not full PR coverage).

## Reproduce
fstar.exe --smt /tmp/gh-aw/fstar/lib/fstar/z3-4.13.3/bin/z3 --z3version 4.13.3 --z3rlimit 10 --cache_dir cache --odir out --include fstar Z3ArithAxiomAudit.fst
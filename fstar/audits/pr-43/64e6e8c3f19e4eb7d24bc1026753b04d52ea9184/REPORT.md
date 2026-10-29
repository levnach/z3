# F* audit PR 43, head 64e6e8c3 (base c95b61cb). Partial scope, semantic audit only.

Proved (Z3Pr43TseitinConstants.fst):
1. tseitin is_marked/is_complement with constants false,(not true),true: sound under every assignment falsifying the clause; new first pass accepts only valid clauses.
2. ite-equality rule (= (ite c u v) y) with equiv y u, equiv y v.
3. euf value/are_equal_args, assuming union-find only relates terms equal in the hint's models (base invariant, not re-proved).

Model: hash-consed terms = structural equality; marks = lists. Standard F* only, no admits.
Excluded: per-gate rules (and/or/eq/ite/implies/xor), xor_holds, main_expr selection, tests/docs, C++ execution/memory safety.
Reproduce: fstar.exe --smt z3-4.13.3/bin/z3 --z3version 4.13.3 --z3rlimit 10 --cache_dir cache --odir out Z3Pr43TseitinConstants.fst (F* 2026.09.27, Z3 4.13.3).
sha256 head tseitin 443871e2..a7707c, head euf 69f183b1..fa846; base tseitin 2b4912b8..5753d, base euf 4ed93662..bd61f.

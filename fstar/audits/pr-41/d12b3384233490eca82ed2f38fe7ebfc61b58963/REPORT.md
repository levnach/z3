# F* audit of PR #41 @ d12b3384233490eca82ed2f38fe7ebfc61b58963 (base 213bdc320dbeb676a240e958d5c0307e7e167b58)

Scope is partial: one changed function, two obligations. Tools: F* 2026.09.27, Z3 4.13.3 (pinned release), z3rlimit 10, fresh cache.

## Audited change
src/sat/sat_proof_trim.cpp, clause-adding routine of the trimmer: the new `std::stable_partition(m_watched, value(lit) != l_false)` before `s.mk_clause`, so the two watched literals are non-false when possible (fixes a missed conflict when both watches are false at base level).

## Obligations (WatchPartition.fst)
1. o1_permutation - PROVED: the partition preserves length and multiset of literals (no literal dropped/duplicated, so the clause added is logically the same).
2. o2_watches - PROVED: if n <= #non-false literals, the first n positions are non-false (n=2: both watches non-false when at least two exist; n=1 otherwise).
(Stability / relative order within groups is not proved; it is not needed for soundness.)

## Correspondence and assumptions
Literals are abstracted to bool (true = value(lit) != l_false). std::stable_partition is modelled by a recursive split + append, not by libstdc++; its standard-guaranteed postcondition is assumed to match. No admit/assume/axioms; only F* standard library (List.Tot) is used.

## Not covered (excluded)
- new throw conditions in trim() (m_has_conflict / theory-axiom empty clause) and the empty-clause handling/set_conflict/propagate(false) logic in the same routine;
- euf_proof_checker has_rule and proof_cmds.cpp guard change (an unsupported-hint hint now treated as an assumed clause: the soundness of assuming such clauses in trimmed proofs is a manual-review item);
- tests, Python examples, docs; C++ execution, memory safety, and the whole solver are not mechanized.

## Reproduce
cd /tmp/gh-aw/agent/fstar-audit && timeout 120 /tmp/gh-aw/fstar/bin/fstar.exe --smt /tmp/gh-aw/fstar/lib/fstar/z3-4.13.3/bin/z3 --z3version 4.13.3 --include fstar --z3rlimit 10 --cache_dir cache --odir out WatchPartition.fst

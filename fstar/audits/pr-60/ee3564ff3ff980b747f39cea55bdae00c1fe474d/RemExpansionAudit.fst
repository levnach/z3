module RemExpansionAudit
(* Model of expand_definitions for (rem u k), k a nonzero integer numeral.
   d = (div u k), m = (mod u k) satisfy the SMT-LIB defining identities. *)
let defs (u k d m:int) : prop = k <> 0 /\ u == k * d + m /\ 0 <= m /\ m < (if k >= 0 then k else -k)
let rem_def (k m:int) : int = if k >= 0 then m else -m
let sign (k:int) : int = if k < 0 then -1 else 1
(* The PR's contribution of c*(rem u k) after expansion: c*s*u - c*s*k*d *)
let expanded (c u k d:int) : int = c * sign k * u - c * sign k * k * d

(* O1: expansion equals c * rem for every coefficient, both signs of k *)
let o1 (c u k d m:int)
  : Lemma (requires defs u k d m) (ensures expanded c u k d == c * rem_def k m) = ()

(* O2: the rem value lies in the range that the accepted clauses imply:
   |rem| < |k|, sign matches k's sign class *)
let o2 (u k d m:int)
  : Lemma (requires defs u k d m)
          (ensures (k > 0 ==> 0 <= rem_def k m /\ rem_def k m < k) /\
                   (k < 0 ==> k < rem_def k m /\ rem_def k m <= 0)) = ()

(* O3: the wrong-sign identity (rem = mod) is not implied for k<0 once m>0:
   concrete witness u=5, k=-3: d=-1, m=2, rem=-2 <> mod=2; the PR rejects it. *)
let o3 () : Lemma (defs 5 (-3) (-1) 2 /\ rem_def (-3) 2 == -2 /\ rem_def (-3) 2 <> 2) = ()

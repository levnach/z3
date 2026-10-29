module Z3ArithAxiomAudit
(* Semantic audit of PR #44, src/sat/smt/arith_theory_checker.h at head 3944722f.
   Values are mathematical integers; this is not a model of C++ execution. *)

let abs_i (n:int) : int = if n < 0 then -n else n

(* SMT-LIB integer div/mod by nonzero k: u = k*d + m, 0 <= m < |k|. *)
let div_mod_spec (u k d m:int) : prop =
  k <> 0 /\ u == k * d + m /\ 0 <= m /\ m < abs_i k

(* O1: expand_definitions: replacing c*(mod u k) by c*u - c*n*(div u k) preserves value. *)
let o1_expand_mod (u n d m c:int)
  : Lemma (requires div_mod_spec u n d m)
          (ensures c * m == c * u - (c * n) * d)
  = ()

(* O1b: replacing c*(u/k) by (c/k)*u is value preserving; scaled by k: k*(u/k)=u. *)
let o1_expand_div (u k q c:int)
  : Lemma (requires k <> 0 /\ k * q == u)
          (ensures k * (c * q) == c * u)
  = ()

(* O2: is_infeasible. Row is c*M + const (<= 0, or < 0 if strict), M = mod u n.
   The checker substitutes M = |n|-1 when c<0 and M=0 otherwise, then check_ineq. *)
let check_ineq (const:int) (strict:bool) : bool =
  const > 0 || (strict && const = 0)

let o2_infeasible (u n d m c const:int) (strict:bool)
  : Lemma (requires div_mod_spec u n d m /\
                    check_ineq (if c < 0 then const + c * (abs_i n - 1) else const) strict)
          (ensures ~ (if strict then c * m + const < 0 else c * m + const <= 0))
  = ()

(* O3: abs case. Literal j's negation is the row r (value rv, <=0 / <0 if strict).
   t = s*u < 0 is added with positive multiplier lam; all coefficients cancel and the
   constant cst is not negative. Then not-L_j and s*u<0 are jointly impossible, and
   abs u = s*u whenever s*u >= 0 (s = 1 or -1), so the clause (abs u = s*u) \/ L_j is valid. *)
let o3_cancel (rv tv lam cst:int) (strict:bool)
  : Lemma (requires lam > 0 /\ cst >= 0 /\ rv + lam * tv == cst /\
                    (if strict then rv < 0 else rv <= 0))
          (ensures ~ (tv < 0))
  = ()

let o3_abs (u s:int)
  : Lemma (requires s == 1 \/ s == -1)
          (ensures s * u >= 0 ==> abs_i u == s * u)
  = ()

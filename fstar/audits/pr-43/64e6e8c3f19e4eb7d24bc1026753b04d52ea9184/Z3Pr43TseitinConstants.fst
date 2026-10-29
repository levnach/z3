module Z3Pr43TseitinConstants
(* Semantic model of PR 43 (levnach/z3, head 64e6e8c3) for
   tseitin::theory_checker::{is_marked,is_complement,check first pass},
   the ite-equality rule, and euf eq_theory_checker::{value,are_equal_args}.
   Terms are hash-consed in C++, so structural equality models pointer equality. *)
open FStar.List.Tot

type bexp =
  | BTrue : bexp
  | BFalse : bexp
  | BVar : nat -> bexp
  | BNot : bexp -> bexp

type iexp =
  | IConst : int -> iexp
  | IVar : nat -> iexp
  | IAdd : iexp -> iexp -> iexp
  | IIte : bexp -> iexp -> iexp -> iexp

let rec beval (s:nat -> bool) (e:bexp) : bool =
  match e with
  | BTrue -> true
  | BFalse -> false
  | BVar n -> s n
  | BNot a -> not (beval s a)

let rec ieval (s:nat -> bool) (t:nat -> int) (e:iexp) : int =
  match e with
  | IConst n -> n
  | IVar n -> t n
  | IAdd a b -> ieval s t a + ieval s t b
  | IIte c a b -> if beval s c then ieval s t a else ieval s t b

(* clause = list of literals; s falsifies the clause *)
let rec falsifies (s:nat -> bool) (l:list bexp) : prop =
  match l with
  | [] -> True
  | x :: r -> beval s x == false /\ falsifies s r

(* complement_mark: mark l, and nmark a when l = (not a) *)
let rec in_mark (l:list bexp) (e:bexp) : bool =
  match l with [] -> false | x :: r -> x = e || in_mark r e
let rec in_nmark (l:list bexp) (e:bexp) : bool =
  match l with
  | [] -> false
  | BNot a :: r -> a = e || in_nmark r e
  | _ :: r -> in_nmark r e

(* tseitin_theory_checker.h *)
let is_false_value (a:bexp) : bool =
  a = BFalse || (match a with BNot BTrue -> true | _ -> false)
let is_marked (l:list bexp) (a:bexp) : bool = in_mark l a || is_false_value a
let is_complement (l:list bexp) (a:bexp) : bool =
  match a with
  | BNot b -> is_marked l b
  | _ -> in_nmark l a || a = BTrue

let rec in_mark_sound (s:nat -> bool) (l:list bexp) (e:bexp)
  : Lemma (requires falsifies s l /\ in_mark l e) (ensures beval s e == false) =
  match l with x :: r -> if x = e then () else in_mark_sound s r e

let rec in_nmark_sound (s:nat -> bool) (l:list bexp) (e:bexp)
  : Lemma (requires falsifies s l /\ in_nmark l e) (ensures beval s e == true) =
  match l with
  | BNot a :: r -> if a = e then () else in_nmark_sound s r e
  | _ :: r -> in_nmark_sound s r e

(* Obligation 1a: marked => false, complement => true, for every falsifying s *)
let marked_sound (s:nat -> bool) (l:list bexp) (a:bexp)
  : Lemma (requires falsifies s l /\ is_marked l a) (ensures beval s a == false) =
  if in_mark l a then in_mark_sound s l a else ()

let complement_sound (s:nat -> bool) (l:list bexp) (a:bexp)
  : Lemma (requires falsifies s l /\ is_complement l a) (ensures beval s a == true) =
  match a with
  | BNot b -> marked_sound s l b
  | _ -> if in_nmark l a then in_nmark_sound s l a else ()

(* Obligation 1b: first pass of check(): a literal that is_complement
   makes the clause valid, i.e. no assignment falsifies it. *)
let rec exists_complement (l0 l:list bexp) : bool =
  match l with [] -> false | x :: r -> is_complement l0 x || exists_complement l0 r

let rec mem_falsified (s:nat -> bool) (m:list bexp) (x:bexp)
  : Lemma (requires falsifies s m /\ mem x m) (ensures beval s x == false) =
  match m with y :: t -> if y = x then () else mem_falsified s t x

let rec first_pass_sound (s:nat -> bool) (l0 l:list bexp)
  : Lemma (requires falsifies s l0 /\ (forall x. mem x l ==> mem x l0) /\ exists_complement l0 l)
          (ensures False) =
  match l with
  | x :: r ->
    if is_complement l0 x then begin
      complement_sound s l0 x;
      assert (mem x l0);
      mem_falsified s l0 x
    end else first_pass_sound s l0 r

let first_pass_valid (l:list bexp)
  : Lemma (requires exists_complement l l) (ensures forall (s:nat -> bool). ~(falsifies s l)) =
  FStar.Classical.forall_intro (fun (s:nat -> bool) ->
    FStar.Classical.move_requires (fun () -> first_pass_sound s l l) ())

(* Obligation 2: new rule (or (= (ite c u v) y) ...) when equiv y u && equiv y v.
   equiv: syntactic, or swapped arguments of the commutative IAdd. *)
let equiv (a b:iexp) : bool =
  a = b || (match a, b with IAdd p q, IAdd p' q' -> p = q' && q = p' | _ -> false)

let equiv_sound (s:nat -> bool) (t:nat -> int) (a b:iexp)
  : Lemma (requires equiv a b) (ensures ieval s t a == ieval s t b) = ()

let ite_eq_sound (s:nat -> bool) (t:nat -> int) (c:bexp) (u v y:iexp)
  : Lemma (requires equiv y u /\ equiv y v)
          (ensures ieval s t (IIte c u v) == ieval s t y) =
  equiv_sound s t y u; equiv_sound s t y v

(* Obligation 3: euf value / are_equal_args.
   eqc models are_equal (union-find). Hypothesis: eqc only relates terms equal
   in every model s of the hint (this is the base invariant of the checker's merges). *)
type lb = | LTrue | LFalse | LUndef

let rec value (eqc:bexp -> bexp -> bool) (e:bexp) (sign:bool) (fuel:nat) : Tot lb (decreases fuel) =
  if eqc e BTrue then (if sign then LFalse else LTrue)
  else if eqc e BFalse then (if sign then LTrue else LFalse)
  else match e with
    | BNot a -> if fuel = 0 then LUndef else value eqc a (not sign) (fuel - 1)
    | _ -> LUndef

let rec size (e:bexp) : nat = match e with BNot a -> 1 + size a | _ -> 0

let are_equal_args (eqc:bexp -> bexp -> bool) (x y:bexp) : bool =
  eqc x y ||
  (let v = value eqc x false (size x) in
   v <> LUndef && v = value eqc y false (size y))

let rec value_sound (eqc:bexp -> bexp -> bool) (s:nat -> bool) (e:bexp) (sign:bool) (fuel:nat)
  : Lemma (requires (forall a b. eqc a b ==> beval s a == beval s b) /\ size e <= fuel)
          (ensures (value eqc e sign fuel == LTrue ==> beval s e == not sign) /\
                   (value eqc e sign fuel == LFalse ==> beval s e == sign))
          (decreases fuel) =
  if eqc e BTrue || eqc e BFalse then ()
  else match e with
    | BNot a -> if fuel = 0 then () else value_sound eqc s a (not sign) (fuel - 1)
    | _ -> ()

let are_equal_args_sound (eqc:bexp -> bexp -> bool) (s:nat -> bool) (x y:bexp)
  : Lemma (requires (forall a b. eqc a b ==> beval s a == beval s b) /\ are_equal_args eqc x y)
          (ensures beval s x == beval s y) =
  if eqc x y then () else begin
    value_sound eqc s x false (size x);
    value_sound eqc s y false (size y)
  end

(* Witness checks: the constants (not true), (not (= true b)) shapes from the PR. *)
let _ = assert_norm (is_marked [] (BNot BTrue) = true)
let _ = assert_norm (is_complement [] BTrue = true)
let _ = assert_norm (is_complement [] (BNot BFalse) = true)
let _ = assert_norm (exists_complement [BVar 0; BNot (BVar 0)] [BVar 0; BNot (BVar 0)] = true)
(* falsifiable look-alike is rejected by the model *)
let _ = assert_norm (exists_complement [BVar 0; BVar 1] [BVar 0; BVar 1] = false)

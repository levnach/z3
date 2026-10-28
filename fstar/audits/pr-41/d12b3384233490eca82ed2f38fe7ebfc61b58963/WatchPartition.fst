module WatchPartition
(* Model of std::stable_partition(m_watched, value(lit) != l_false) in
   sat_proof_trim.cpp. A literal is abstracted to its predicate value:
   true = not false under the base-level assignment. *)
open FStar.List.Tot

let rec split (l:list bool) : list bool & list bool =
  match l with
  | [] -> ([], [])
  | x::t -> let (a,b) = split t in if x then (x::a, b) else (a, x::b)

let partition (l:list bool) : list bool = let (a,b) = split l in a @ b

let rec all_eq (v:bool) (l:list bool) : bool =
  match l with [] -> true | x::t -> (x = v) && all_eq v t

let rec split_props (l:list bool)
  : Lemma (let (a,b) = split l in
           all_eq true a /\ all_eq false b /\
           length a == count true l /\ length b == count false l)
  = match l with [] -> () | _::t -> split_props t

let rec cnt (v:bool) (x y:list bool)
  : Lemma (count v (x @ y) == count v x + count v y) =
  match x with [] -> () | _::t -> cnt v t y

let rec ac (v:bool) (x:list bool)
  : Lemma (requires all_eq v x) (ensures count v x == length x /\ count (not v) x == 0) =
  match x with [] -> () | _::t -> ac v t

let rec count_sum (l:list bool)
  : Lemma (count true l + count false l == length l) =
  match l with [] -> () | _::t -> count_sum t

(* O1: the partition is a permutation: same length and same multiset of truth values. *)
let o1_permutation (l:list bool)
  : Lemma (let r = partition l in
           length r == length l /\ count true r == count true l /\ count false r == count false l)
  = split_props l; count_sum l;
    let (a,b) = split l in
    FStar.List.Tot.Properties.append_length a b;
    cnt true a b; cnt false a b;
    ac true a; ac false b

(* O2: the first n positions hold non-false literals whenever n <= #non-false,
   so both watched literals are non-false if at least two exist. *)
let rec prefix_true (n:nat) (l:list bool) : bool =
  match n, l with
  | 0, _ -> true
  | _, [] -> false
  | _, x::t -> x && prefix_true (n-1) t

let rec prefix_append (n:nat) (a b:list bool)
  : Lemma (requires all_eq true a /\ n <= length a)
          (ensures prefix_true n (a @ b))
  = match n, a with 0, _ -> () | _, _::t -> prefix_append (n-1) t b

let o2_watches (l:list bool) (n:nat)
  : Lemma (requires n <= count true l)
          (ensures prefix_true n (partition l))
  = split_props l;
    let (a,b) = split l in
    prefix_append n a b

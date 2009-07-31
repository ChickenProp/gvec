(ns net.philh.gvec
  (:refer-clojure :exclude (+ - * /))
  (:use clojure.contrib.generic
	clojure.contrib.generic.arithmetic
	clojure.contrib.types
	clojure.contrib.import-static))

(import-static java.lang.Math PI sqrt sin cos atan atan2)

(declare length)
(declare extend)
(declare with-dim)
(declare vec-to-gvec)
(declare keyword-to-index)

;; If someone doesn't have java.awt.Point, this will error. Which is a shame,
;; because there's no reason to require them to have it. Best way to avoid that?
(deftype ::gvec gvec
  "Returns a gvec (geometrical vector) with the appropriate indices.

Indices can be specified in the following ways:
Cartesian: (gvec i_1 ... i_n), (gvec [i_1 ... i_n])
  In the second form, any seqable argument works except a map.
Polar: (gvec :r len :theta ang :phi inc), (gvec {:r len :theta ang :phi inc})
  :r defaults to 1, :theta to 0. If :phi is not provided, the gvec returned
  will only have two dimensions, which is almost the same as saying :phi
  defaults to pi/2.  :theta is the azimuth and :phi the incidence, which is
  somewhat unconventional in spherical coordinates but is consistent with 2d
  polar and cylindrical coordinates.
From a java.awt.Point: (gvec point)
Empty: (gvec), (gvec [])
  The second form can take any seqable argument except a map, as long as it's
  empty. I'm not sure why you'd want an empty gvec anyway.
There is currently no support for generalised n-dimensional hyperspherical
coordinates."
  (fn
    ([] (vec-to-gvec []))
    ([v]
       (cond (map? v)
	       (let [r (get v :r 1)
		     the (get v :theta 0)
		     phi (get v :phi nil)
		     x (* r (cos the))
		     y (* r (sin the))]
		 (vec-to-gvec (if phi
				[(* x (sin phi))
				 (* y (sin phi))
				 (* r (cos phi))]
				[x y])))
	     (instance? java.awt.Point v)
	       (vec-to-gvec [(.getX #^java.awt.Point v)
			     (.getY #^java.awt.Point v)])
	     (instance? clojure.lang.Seqable v) (vec-to-gvec v)
	     true (vec-to-gvec [v])))
    ([v1 & v+] ; it would be nice if we could use :as here
       (let [vs (cons v1 v+)]
	 (if (keyword? v1)
	   (gvec (apply hash-map vs))
	   (vec-to-gvec vs))))))

(derive ::gvec root-type)

(defn vec-to-gvec
  "Returns a geometrical version of a vector. Does not attach extra metadata."
  [vecable]
  (let [v (vec vecable)]
    (proxy [clojure.lang.APersistentVector] [(meta vecable)]
      (seq [] (seq v))
      (nth [idx] (nth v idx 0))
      (assoc [idx val]
	(let [i (keyword-to-index idx)]
	  (if (> i (count v))	     ;assoc for vectors requires index <= count.
	    (assoc (extend this (+ 1 i)) i val)
	    (gvec (assoc v i val)))))
      (count [] (count v))
      (cons [x] (gvec (conj v x)))
      (invoke [idx]
	(cond (number? idx) (nth v idx 0)
	      (vector? idx)
	        (if (= (idx 0) :pol)
		  (let [idx (idx 1)]
		    (cond (= idx :r) (length v)
			  (= idx 0) (atan2 (nth v 1 0) (nth v 0 0))
			  (= (nth v (inc idx) 0) 0) (/ PI 2)
			  true (atan (/ (length (with-dim this (inc idx)))
					     (nth v (inc idx) 0))))))
	      (keyword? idx) (this (keyword-to-index idx))))
      (withMeta [m] (vec-to-gvec (with-meta v m))))))

(defn gvec-polar
  "Returns a gvec defined in hyperspherical coordinates. Currently works only up
to four dimensions: (r, theta, phi, psi)."
  ([r]
     (gvec r))
  ([r theta]
     (* r (gvec (cos theta) (sin theta))))
  ([r theta phi]
     (* r (gvec (* (cos theta) (sin phi))
		(* (sin theta) (sin phi))
		(cos phi))))
  ([r theta phi psi]
     (* r (gvec (* (cos theta) (sin phi) (sin psi))
		(* (sin theta) (sin phi) (sin psi))
		(* (cos phi) (sin psi))
		(cos psi)))))

(defn keyword-to-index
  "Given a keyword, translates it to the appropriate numeric or vector index
for use with gvecs."
  [kw]
  (let [trans {:x 0, :y 1, :z 2, :w 3,
	       :r [:pol :r], :theta [:pol 0], :phi [:pol 1]}]
    (or (trans kw) kw)))

(defn extend
  "Returns a gvec with the same coordinates as its argument, but with
 (max n (count gv)) dimensions."
  [gv n]
  (gvec (concat (seq gv) (repeat (- n (count gv)) 0))))

(defn with-dim
  "Returns a gvec with the same coordinates as its argument, but with n
dimensions. If the original has more than n dimensions, the excess ones will be
lost."
  [gv n]
  (cond (> n (count gv))
	  (extend gv n)
	(< n (count gv))
	  (gvec (take n gv))
	true gv))

(defn gvmap
  "Maps over one or more gvecs, returning another gvec. Normal map would stop
early if one of the arguments had fewer dimensions than the others."
  [f & gvs]
  (let [max-count (apply max (map count gvs))]
    (gvec (apply map f
		 (map #(extend % max-count)
		      gvs)))))

(defmethod + [::gvec ::gvec]
  [a b]
  (gvmap + a b))
(defmethod - ::gvec
  [a]
  (gvmap - a))
(defmethod * [::gvec Number]
  [a n]
  (gvmap #(* % n) a))
(defmethod * [Number ::gvec]
  [n a]
  (* a n))

(defn length
  [a]
  (sqrt (apply + (map #(* % %) a))))

(defn normalise
  [a]
  (let [l (length a)]
    (if (= l 0) a
	(/ a (length a)))))

(defn plane-rotate
  "Rotates a point anticlockwise about the origin on the ax1-ax2 plane."
  [gv axis-1 axis-2 angle]
  (let [g (gvec gv)
	v1 (g axis-1)
	v2 (g axis-2)
	c (cos angle)
	s (sin angle)]
    (assoc g
      axis-1 (- (* v1 c) (* v2 s))
      axis-2 (+ (* v1 s) (* v2 c)))))


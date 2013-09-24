(ns folding.core
  (:gen-class))

(def propensities {"V" [0.91, 2.00],
                   "I" [1.04, 1.79],
                   "L" [1.28, 1.15],
                   "M" [1.26, 1.01],
                   "P" [0.44, 0.40],
                   "A" [1.41 0.75],
                   "C" [0.85 1.36],
                   "F" [1.00 1.40],
                   "Y" [0.98 1.37],
                   "W" [1.07 1.23],
                   "Q" [1.26 0.72],
                   "S" [0.76 0.81],
                   "T" [0.78 1.21],
                   "N" [0.73, 0.63],
                   "H" [0.87 0.99],
                   "D" [0.82 0.55],
                   "K" [1.17 0.76],
                   "E" [1.39 0.65],
                   "R" [1.21 0.85],
                   "G" [0.44 0.67]})

(defn not-in?
  "returns true if a collection does not conatin a specific entry"
  [entry collection]
  (nil? (some #(= entry %) collection)))

(defn has?
  "returns true if a collection contains a specific entry"
  [entry collection]
  (not (not-in? entry collection)))

(defn char-at
  "character at position"
  [position string]
  (str (get string position)))

(defn get-propensity
  "gets the propensity for an amino acid to paticipate in an alpha-helix or a beta-sheet"
  [character-position protein-string]
  (get propensities (char-at character-position protein-string)))

(defn get-alpha-propensity
  "gets the propensity for an amino-acid to participate in an alpha helix"
  [character-position protein-string]
  (get (get-propensity character-position protein-string) 0))

(defn get-beta-propensity
  "gets the propensity for an amino acid to participate in an alpha-helix"
  [character-position protein-string]
  (get (get-propensity character-position protein-string) 1))

(defn has-propensity-for
  "given a propensity function, calculates if two amino acids would bind"
  [propensity-function row column protein-string]
  (>= (* (propensity-function row protein-string) (propensity-function column protein-string))  1.1))

(defn amino-acid-propensities
  "Calculates the binding of an amino acid to other amino acids in the protein string. This function returns a vector of all contact points that favour a binding."
  [propensity-function row protein-string]
  (let [length (dec (count protein-string))]
    (vec (for [column (range row length) :when (has-propensity-for propensity-function row column protein-string)] [row column])) ))

(defn propensities-for
  "calculates a propensity row"
  [row protein-string]
  [(amino-acid-propensities get-alpha-propensity row protein-string) (amino-acid-propensities get-beta-propensity row protein-string)])

(defn initial-map
  "creates initial map"
  [protein-string]
  (loop [i 0
         alpha-helixes []
         beta-sheets []]
    (if (> i (dec (count protein-string)))
      [alpha-helixes beta-sheets]
      (let [calculated-propensities (propensities-for i protein-string)
            alpha-propensities (get calculated-propensities 0)
            beta-propensities (get calculated-propensities 1)]
        (recur (inc i) (concat alpha-helixes alpha-propensities) (concat beta-sheets beta-propensities))
        ))))

(defn possible-alpha?
  "could this amino-acid combination be part of an alpha helix"
  [[i j] protein-string]
  (let [row-propensity (get-alpha-propensity i protein-string)
        column-propensity (get-alpha-propensity j protein-string)]
    (or (>= row-propensity 1.1)
        (>= column-propensity 1.1)
        (and (>= row-propensity 0.9)
             (>= column-propensity 0.9))))
  )

(defn get-alpha-neighbours
  "returns the nearest neighbours on a alph-helix"
  [[i j] length]
  (filter (fn [[i j]] (and (>= i 0) (>= j 0) (<= i length) (<= j length))) [[(dec i) (dec j)] [(inc i) (inc j)]])
  )

(defn could-be-alpha
  ""
  [known-alpha protein-string]
  (let [protein-length (count protein-string)]
    (filter (fn [neighbour] (possible-alpha? neighbour protein-string)) (get-alpha-neighbours known-alpha protein-length)))
  )

(defn print-join
  "prints a join"
  [type [i j]]
  (println (str i " " j " " type)))

(defn explore-for-more-alphas
  "tries to build an alpha helix"
  [definite-alphas protein-string]
  (loop [old definite-alphas]
    (let [discovered (reduce concat (map (fn [known-alpha] (could-be-alpha known-alpha protein-string)) definite-alphas))
          count-old (count old)
          combined (set (concat old discovered))
          count-set (count combined)]
      (if (= count-set count-old)
        old
        (recur combined)))))

(defn in-same-parallel-grid?
  "checks if entries are near"
  [distance [i j] [k l]]
  (and (>= k i)
       (>= l j)
       (> (+ i distance) k)
       (> (+ j distance) l))
  )

(defn in-same-anti-parallel-grid?
  "checks if entries are near"
  [distance [i j] [k l]]
  (and (>= k i)
       (>= j l)
       (> (+ i distance) k)
       (> (+ l distance) j))
  )

(defn exists-in?
  "checks for existence at relative position"
  [betas candidate]
  (has? candidate betas))

(defn parallel-template
  "a template for parallel beta sheets"
  [[i j]]
  [[i j] [(+ i 2) j] [(+ i 2) (+ j 2)] [(+ i 4) (+ j 2)] [(+ i 4) (+ j 4)]]
  )

(defn valid-parallel?
  "can be a parallel arrangement"
  [betas]
  (let [ exists-in-betas? (partial exists-in? betas)]
    (reduce #(and %1 %2) (map exists-in-betas? (parallel-template (first betas)))))
  )

(defn create-parallel-from
  "removes joins that would make the matrix invalid as a parallel beta sheet"
  [betas]
  (parallel-template (first betas)))

(defn anti-parallel-template
  "a template for anti-parallel beta sheets"
  [[i j]]
  [[i j] [(+ i 2) (- j 2)] [(+ i  4) (- j 4)]])

(defn valid-anti-parallel?
  "can be an anti-parallel arrangement"
  [betas]
  (let [exists-in-betas? (partial exists-in? betas)]
    (reduce #(and %1 %2) (map exists-in-betas? (anti-parallel-template (first betas))))
    )
  )

(defn all-potential-beta-grids
  "Uses the provided function to group all potential participants in a beta sheet for a 5x5 matrix"
  [betas grid-fn]
  (map (fn [beta] (filter (partial grid-fn 5 beta) betas)) betas)
  )

(defn in-a-sequence?
  [[x y] alphas]
  (let [in-sequence? (fn [alpha] (has? alpha alphas))]
    (or (and (in-sequence? [(inc x) (inc y)])
             (in-sequence? [(+ x 2) (+ y 2)])
             (in-sequence? [(+ x 3) (+ y 3)]))
        (and (in-sequence? [(dec x) (dec y)])
             (in-sequence? [(inc x) (inc y)])
             (in-sequence? [(+ x 2) (+ y 2)]))
        (and (in-sequence? [(- x 2) (- y 2)])
             (in-sequence? [(dec x) (dec y)])
             (in-sequence? [(+ x 1) (+ y 1)]))
        (and (in-sequence? [(- x 3) (- y 3)])
             (in-sequence? [(- x 2) (- y 2)])
             (in-sequence? [(dec x) (dec y)])))))


(defn strip-short-sequences
  "removes all alphas that are not part of a sequence at least 4 amino acids long"
  [alphas]
  (filter (fn [alpha] (in-a-sequence? alpha alphas)) alphas))

(defn can-exist?
  "can the contact exist given other known contacts"
  [[x y] entries]
  (empty? (for [[i j] entries
                :when (or (= x i)
                          (= y j))] [i j])))

(defn create-map
  "creates a contact map by getting an initial map based on propensity and then refining it"
  [protein-string]
  ;; create the initial lists of alphas and betas based on propensities
  (let [[first-alphas first-betas] (initial-map protein-string)
        ;; ensure that alphas have a 4 amino-acid offset
        definite-alphas (filter (fn [[x y]] (= 4 (- y x))) first-alphas)
        ;; given the definite alphas, explore the neighbourhood for ones that could participate
        possible-alphas (set (explore-for-more-alphas definite-alphas protein-string))
        ;; remove any sequences of alphas that are less than 4 amino-acids long
        alphas (strip-short-sequences possible-alphas)
        ;; remove any illegal betas (a gap of less than 3 amino acids)
        possible-betas (filter (fn [[x y]] (> (- y x) 2)) first-betas)
        ;; remove any betas that fall on the same x or y coordinates as the alphas
        betas-ruled-out-by-alphas (filter (fn [beta] (can-exist? beta alphas)) possible-betas)
        ;; create all groupings of all betas in 5x5 grids anchored by a known beta in the top left corner (for parallel beta sheets)
        all-potential-parallel-grids (all-potential-beta-grids betas-ruled-out-by-alphas in-same-parallel-grid?)
        ;; create all groupings of all betas in 5x5 grids anchored by a known beta in the top right hand corner (for parallel beta sheets)
        all-potential-anti-parallel-grids (all-potential-beta-grids betas-ruled-out-by-alphas in-same-anti-parallel-grid?)
        ;; for each parallel grouping, eliminate the betas that do not conform to a parallel beta sheet
        refined-parallel-betas (set (reduce concat (map create-parallel-from (filter valid-parallel? all-potential-parallel-grids))))
        ;; for each anti-parallel grouping, eliminate the betas that do not conform to an anti-parallel beta sheet
        refined-anti-parallel-betas (set (reduce concat (filter valid-anti-parallel? all-potential-anti-parallel-grids)))
        ] [alphas (vec (concat refined-parallel-betas refined-anti-parallel-betas))]))


(defn contact-map-print
  "creates a string for printing the contact map"
  [[alphas betas]]
  (let [print-a (partial print-join "A")
        print-b (partial print-join "B")]
    (map print-a alphas)
    (map print-b betas)))

(defn -main
  "Entry point to the contact-map generation program."
  [& args]
  ;; work around dangerous default behaviour in Clojure
  (alter-var-root #'*read-eval* (constantly false))
  (let [
        stripped (clojure.string/replace (first args) #"\s" "")
        [alphas betas] (create-map stripped)
        ]
    (println (count stripped))
    (doseq [[i j] alphas] (println (str i " " j " A")))
    (doseq [[i j] betas] (println (str i " " j " B")))
    )
  )

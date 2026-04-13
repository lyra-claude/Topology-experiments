module Main where

import IslandGA (Topology, Stats(..), runSimulation)
import Maze (Maze)    -- Domain instance for 15x15 Maze
import Maze8 (Maze8)  -- Domain instance for 8x8 Maze
import OneMax (OneMax) -- Domain instance for OneMax
import NKLandscape (NK0Individual, NK2Individual, NK4Individual, NK6Individual)
import Domain (Domain(..))

import qualified Data.Vector as V
import Data.Proxy (Proxy(..))
import System.Environment (getArgs)
import System.Random (mkStdGen)
import Data.Bits (xor, testBit, popCount)
import System.IO (hFlush, stdout, hPutStrLn, stderr)

-- ---------------------------------------------------------------------------
-- Topology builders
-- ---------------------------------------------------------------------------

-- | No edges -- fully disconnected islands.
disconnected :: Int -> Topology
disconnected n = V.replicate n []

-- | Ring (cycle): each island connected to its two neighbors.
ring :: Int -> Topology
ring n = V.generate n (\i -> [(i - 1) `mod` n, (i + 1) `mod` n])

-- | Star: island 0 connected to all others, others connected only to 0.
star :: Int -> Topology
star n = V.generate n (\i -> if i == 0
                               then [1 .. n - 1]
                               else [0])

-- | Complete graph: every island connected to every other.
complete :: Int -> Topology
complete n = V.generate n (\i -> [j | j <- [0 .. n - 1], j /= i])

-- | Hypercube of dimension k (n = 2^k nodes).
-- Two nodes are connected iff they differ in exactly one bit.
hypercube :: Int -> Topology
hypercube k =
  let n = 2 ^ k
  in V.generate n (\i -> [j | j <- [0 .. n - 1], popCount (i `xor` j) == 1])

-- | Barbell: two cliques of n/2 connected by a single edge.
-- Nodes 0..half-1 form clique 1, nodes half..n-1 form clique 2.
-- Node (half-1) and node half are the bridge.
barbell :: Int -> Topology
barbell n =
  let half = n `div` 2
  in V.generate n (\i ->
       if i < half
         then
           -- Clique 1: connected to all other clique-1 nodes
           let clique = [j | j <- [0 .. half - 1], j /= i]
           in if i == half - 1
                then half : clique  -- bridge node
                else clique
         else
           -- Clique 2: connected to all other clique-2 nodes
           let clique = [j | j <- [half .. n - 1], j /= i]
           in if i == half
                then (half - 1) : clique  -- bridge node
                else clique
     )

-- | Watts-Strogatz small-world graph.
wattsStrogatz :: Int -> Int -> Double -> Int -> Topology
wattsStrogatz n k p seed =
  let halfK = k `div` 2
      ringLattice = V.generate n (\i ->
        [j | d <- [1..halfK], let j = (i + d) `mod` n] ++
        [j | d <- [1..halfK], let j = (i - d) `mod` n])

      rewire :: V.Vector [Int] -> V.Vector [Int]
      rewire adj0 = foldl rewireEdge adj0
        [(i, (i + d) `mod` n) | i <- [0..n-1], d <- [1..halfK]]
        where
          shouldRewire i j =
            let h = (i * 7919 + j * 6271 + seed * 1031) `mod` 10000
            in fromIntegral h / 10000.0 < p

          rewireEdge :: V.Vector [Int] -> (Int, Int) -> V.Vector [Int]
          rewireEdge adj (i, j)
            | not (shouldRewire i j) = adj
            | otherwise =
                let newJ = ((i * 3571 + j * 2749 + seed * 947) `mod` (n - 1))
                    newJ' = if newJ >= i then newJ + 1 else newJ
                in if newJ' == j || newJ' `elem` (adj V.! i)
                     then adj
                     else
                       let adjI = filter (/= j) (adj V.! i) ++ [newJ']
                           adjJ = filter (/= i) (adj V.! j)
                           adjN = (adj V.! newJ') ++ [i]
                           adj' = adj V.// [(i, adjI), (j, adjJ), (newJ', adjN)]
                       in adj'

  in rewire ringLattice

-- | Random 3-regular graph (deterministic construction).
randomRegular :: Int -> Int -> Int -> Topology
randomRegular n d seed =
  let base = ring n
      addExtras :: V.Vector [Int] -> V.Vector [Int]
      addExtras adj = foldl tryAddEdge adj [0 .. n - 1]
        where
          tryAddEdge :: V.Vector [Int] -> Int -> V.Vector [Int]
          tryAddEdge adj' i =
            let currentDeg = length (adj' V.! i)
            in if currentDeg >= d
                 then adj'
                 else
                   let target = ((i * 5021 + seed * 1733) `mod` (n - 2))
                       target' = if target >= i then target + 1 else target
                       targetDeg = length (adj' V.! target')
                   in if target' `elem` (adj' V.! i) || targetDeg >= d
                        then adj'
                        else adj' V.// [ (i, target' : (adj' V.! i))
                                       , (target', i : (adj' V.! target'))
                                       ]
  in addExtras base

-- ---------------------------------------------------------------------------
-- Bridge experiment: iso-spectral families with varying beta_1
-- ---------------------------------------------------------------------------

-- | Add an undirected edge (u,v) to an adjacency-list topology.
addEdge :: Int -> Int -> Topology -> Topology
addEdge u v topo =
  topo V.// [ (u, v : (topo V.! u))
            , (v, u : (topo V.! v))
            ]

-- Family 1: lambda_2 = 0.5858 (cycle-based, chords between equal-Fiedler vertices)

-- | Ring + 1 chord: C_8 with edge (0,4). beta_1=2, lambda_2=0.5858.
ringChord1 :: Topology
ringChord1 = addEdge 0 4 (ring 8)

-- | Ring + 2 chords: C_8 with edges (0,4),(1,3). beta_1=3, lambda_2=0.5858.
ringChord2 :: Topology
ringChord2 = addEdge 1 3 (addEdge 0 4 (ring 8))

-- | Ring + 3 chords: C_8 with edges (0,4),(1,3),(5,7). beta_1=4, lambda_2=0.5858.
ringChord3 :: Topology
ringChord3 = addEdge 5 7 (addEdge 1 3 (addEdge 0 4 (ring 8)))

-- Family 2: lambda_2 = 1.0 (star-based, leaf-leaf edges)

-- | Star + 1 leaf edge: S_8 with edge (1,2). beta_1=1, lambda_2=1.0.
starLeaf1 :: Topology
starLeaf1 = addEdge 1 2 (star 8)

-- | Star + 2 leaf edges: S_8 with edges (1,2),(3,4). beta_1=2, lambda_2=1.0.
starLeaf2 :: Topology
starLeaf2 = addEdge 3 4 (addEdge 1 2 (star 8))

-- | Star + 3 leaf edges: S_8 with edges (1,2),(3,4),(5,6). beta_1=3, lambda_2=1.0.
starLeaf3 :: Topology
starLeaf3 = addEdge 5 6 (addEdge 3 4 (addEdge 1 2 (star 8)))

-- ---------------------------------------------------------------------------
-- Cycle-length experiment: two-cycle-bridge graphs, all beta_1=2
-- Two cycles of length n1, n2 sharing a single bridge edge (node 0 -- node n1).
-- ---------------------------------------------------------------------------

-- | Build a two-cycle bridge graph: cycle of n1 nodes (0..n1-1),
-- cycle of n2 nodes (n1..n1+n2-1), bridge edge between 0 and n1.
twoCycleBridge :: Int -> Int -> Topology
twoCycleBridge n1 n2 =
  let totalN = n1 + n2
      -- Cycle 1: nodes 0..n1-1
      cycle1 i = [(i - 1) `mod` n1, (i + 1) `mod` n1]
      -- Cycle 2: nodes n1..n1+n2-1
      cycle2 i =
        let local = i - n1
        in [n1 + ((local - 1) `mod` n2), n1 + ((local + 1) `mod` n2)]
      base = V.generate totalN (\i ->
        if i < n1 then cycle1 i else cycle2 i)
      -- Add bridge edge 0 <-> n1
  in addEdge 0 n1 base

-- | G1: two 3-cycles with bridge. beta_1=2, mean_cycle=3.0, n=6, lambda_2=0.4384
cycleBridge33 :: Topology
cycleBridge33 = twoCycleBridge 3 3

-- | G2: 3-cycle + 9-cycle with bridge. beta_1=2, mean_cycle=6.0, n=12, lambda_2=0.1907
cycleBridge39 :: Topology
cycleBridge39 = twoCycleBridge 3 9

-- | G3: two 9-cycles with bridge. beta_1=2, mean_cycle=9.0, n=18, lambda_2=0.0822
cycleBridge99 :: Topology
cycleBridge99 = twoCycleBridge 9 9

-- ---------------------------------------------------------------------------
-- Directed topology builders (n=8, m=16 directed edges each)
-- Used for density-cycle confound experiment: constant density, varying
-- simple directed cycle count.
-- Adjacency list interpretation: topo V.! j = in-neighbors of j
-- (i.e., islands that send migrants TO island j).
-- ---------------------------------------------------------------------------

-- | DAG-Layer: layered {0,1} -> {2,3,4} -> {5,6,7}, plus 0->1.
-- 0 directed cycles.
dagLayer :: Topology
dagLayer = V.fromList [ []
                      , [0]
                      , [0, 1]
                      , [0, 1]
                      , [0, 1]
                      , [2, 3, 4]
                      , [2, 3, 4]
                      , [2, 3, 4]
                      ]

-- | DAG-Wide: node 0 broadcasts to all, layered fan-out.
-- 0 directed cycles.
dagWide :: Topology
dagWide = V.fromList [ []
                     , [0]
                     , [0]
                     , [0]
                     , [0]
                     , [0, 1, 2, 3]
                     , [0, 1, 2, 3, 4]
                     , [0, 4, 5]
                     ]

-- | LowCyc-1: DAG-Layer with one feedback edge (5->0).
-- Creates cycle 0->2->5->0. 3 directed cycles.
lowcyc1 :: Topology
lowcyc1 = V.fromList [ [5]
                      , []
                      , [0, 1]
                      , [0, 1]
                      , [0, 1]
                      , [2, 3, 4]
                      , [2, 3, 4]
                      , [2, 3, 4]
                      ]

-- | Bidirectional ring: forward and backward directed rings.
-- 10 directed cycles.
bidirRing :: Topology
bidirRing = V.fromList [ [1, 7]
                       , [0, 2]
                       , [1, 3]
                       , [2, 4]
                       , [3, 5]
                       , [4, 6]
                       , [5, 7]
                       , [0, 6]
                       ]

-- | Two disconnected directed 4-cliques (ring + cross-chords each).
-- 14 directed cycles.
twoCliques :: Topology
twoCliques = V.fromList [ [2, 3]
                         , [0, 3]
                         , [0, 1]
                         , [1, 2]
                         , [6, 7]
                         , [4, 7]
                         , [4, 5]
                         , [5, 6]
                         ]

-- | 2x4 mesh with vertical bidirectional edges and horizontal wrap-around.
-- 20 directed cycles.
meshCyclic :: Topology
meshCyclic = V.fromList [ [3, 4]
                        , [0, 5]
                        , [1, 6]
                        , [2, 7]
                        , [0, 7]
                        , [1, 4]
                        , [2, 5]
                        , [3, 6]
                        ]

-- | Overlapping directed triangles: 0->1->2->0, 2->3->4->2, 4->5->6->4,
-- 6->7->0->6, plus skip chain 1->3->5->7->1.
-- 29 directed cycles.
denseTriangles :: Topology
denseTriangles = V.fromList [ [2, 7]
                            , [0, 7]
                            , [1, 4]
                            , [1, 2]
                            , [3, 6]
                            , [3, 4]
                            , [0, 5]
                            , [5, 6]
                            ]

-- | Ring + skip-2 ring: i->i+1 and i->i+2 (mod 8).
-- 47 directed cycles (maximum for n=8, m=16).
ringSkip2 :: Topology
ringSkip2 = V.fromList [ [6, 7]
                        , [0, 7]
                        , [0, 1]
                        , [1, 2]
                        , [2, 3]
                        , [3, 4]
                        , [4, 5]
                        , [5, 6]
                        ]

-- ---------------------------------------------------------------------------
-- Foster census cubic symmetric graphs (13 graphs, all 3-regular)
-- ---------------------------------------------------------------------------

-- | K4 (4 vertices, cubic symmetric)
fosterK4 :: Topology
fosterK4 = V.fromList
  [ [1, 2, 3]
  , [0, 2, 3]
  , [0, 1, 3]
  , [0, 1, 2]
  ]

-- | K33 (6 vertices, cubic symmetric)
fosterK33 :: Topology
fosterK33 = V.fromList
  [ [3, 4, 5]
  , [3, 4, 5]
  , [3, 4, 5]
  , [0, 1, 2]
  , [0, 1, 2]
  , [0, 1, 2]
  ]

-- | Cube / Q3 (8 vertices, cubic symmetric)
fosterCube :: Topology
fosterCube = V.fromList
  [ [1, 3, 4]
  , [0, 2, 7]
  , [1, 3, 6]
  , [0, 2, 5]
  , [0, 5, 7]
  , [3, 4, 6]
  , [2, 5, 7]
  , [1, 4, 6]
  ]

-- | Petersen (10 vertices, cubic symmetric)
fosterPetersen :: Topology
fosterPetersen = V.fromList
  [ [1, 4, 5]
  , [0, 2, 6]
  , [1, 3, 7]
  , [2, 4, 8]
  , [0, 3, 9]
  , [0, 7, 8]
  , [1, 8, 9]
  , [2, 5, 9]
  , [3, 5, 6]
  , [4, 6, 7]
  ]

-- | Heawood (14 vertices, cubic symmetric)
fosterHeawood :: Topology
fosterHeawood = V.fromList
  [ [1, 5, 13]
  , [0, 2, 10]
  , [1, 3, 7]
  , [2, 4, 12]
  , [3, 5, 9]
  , [0, 4, 6]
  , [5, 7, 11]
  , [2, 6, 8]
  , [7, 9, 13]
  , [4, 8, 10]
  , [1, 9, 11]
  , [6, 10, 12]
  , [3, 11, 13]
  , [0, 8, 12]
  ]

-- | Mobius-Kantor (16 vertices, cubic symmetric)
fosterMobiusKantor :: Topology
fosterMobiusKantor = V.fromList
  [ [1, 5, 15]
  , [0, 2, 12]
  , [1, 3, 7]
  , [2, 4, 14]
  , [3, 5, 9]
  , [0, 4, 6]
  , [5, 7, 11]
  , [2, 6, 8]
  , [7, 9, 13]
  , [4, 8, 10]
  , [9, 11, 15]
  , [6, 10, 12]
  , [1, 11, 13]
  , [8, 12, 14]
  , [3, 13, 15]
  , [0, 10, 14]
  ]

-- | Pappus (18 vertices, cubic symmetric)
fosterPappus :: Topology
fosterPappus = V.fromList
  [ [1, 5, 17]
  , [0, 2, 8]
  , [1, 3, 13]
  , [2, 4, 10]
  , [3, 5, 15]
  , [0, 4, 6]
  , [5, 7, 11]
  , [6, 8, 14]
  , [1, 7, 9]
  , [8, 10, 16]
  , [3, 9, 11]
  , [6, 10, 12]
  , [11, 13, 17]
  , [2, 12, 14]
  , [7, 13, 15]
  , [4, 14, 16]
  , [9, 15, 17]
  , [0, 12, 16]
  ]

-- | Dodecahedron (20 vertices, cubic symmetric)
fosterDodecahedron :: Topology
fosterDodecahedron = V.fromList
  [ [1, 10, 19]
  , [0, 2, 8]
  , [1, 3, 6]
  , [2, 4, 19]
  , [3, 5, 17]
  , [4, 6, 15]
  , [2, 5, 7]
  , [6, 8, 14]
  , [1, 7, 9]
  , [8, 10, 13]
  , [0, 9, 11]
  , [10, 12, 18]
  , [11, 13, 16]
  , [9, 12, 14]
  , [7, 13, 15]
  , [5, 14, 16]
  , [12, 15, 17]
  , [4, 16, 18]
  , [11, 17, 19]
  , [0, 3, 18]
  ]

-- | Desargues (20 vertices, cubic symmetric)
fosterDesargues :: Topology
fosterDesargues = V.fromList
  [ [1, 5, 19]
  , [0, 2, 16]
  , [1, 3, 11]
  , [2, 4, 14]
  , [3, 5, 9]
  , [0, 4, 6]
  , [5, 7, 15]
  , [6, 8, 18]
  , [7, 9, 13]
  , [4, 8, 10]
  , [9, 11, 19]
  , [2, 10, 12]
  , [11, 13, 17]
  , [8, 12, 14]
  , [3, 13, 15]
  , [6, 14, 16]
  , [1, 15, 17]
  , [12, 16, 18]
  , [7, 17, 19]
  , [0, 10, 18]
  ]

-- | Nauru / GP(12,5) (24 vertices, cubic symmetric)
fosterNauru :: Topology
fosterNauru = V.fromList
  [ [1, 11, 12]
  , [0, 2, 13]
  , [1, 3, 14]
  , [2, 4, 15]
  , [3, 5, 16]
  , [4, 6, 17]
  , [5, 7, 18]
  , [6, 8, 19]
  , [7, 9, 20]
  , [8, 10, 21]
  , [9, 11, 22]
  , [0, 10, 23]
  , [0, 17, 19]
  , [1, 18, 20]
  , [2, 19, 21]
  , [3, 20, 22]
  , [4, 21, 23]
  , [5, 12, 22]
  , [6, 13, 23]
  , [7, 12, 14]
  , [8, 13, 15]
  , [9, 14, 16]
  , [10, 15, 17]
  , [11, 16, 18]
  ]

-- | F26A (26 vertices, cubic symmetric)
fosterF26A :: Topology
fosterF26A = V.fromList
  [ [1, 19, 25]
  , [0, 2, 8]
  , [1, 3, 21]
  , [2, 4, 10]
  , [3, 5, 23]
  , [4, 6, 12]
  , [5, 7, 25]
  , [6, 8, 14]
  , [1, 7, 9]
  , [8, 10, 16]
  , [3, 9, 11]
  , [10, 12, 18]
  , [5, 11, 13]
  , [12, 14, 20]
  , [7, 13, 15]
  , [14, 16, 22]
  , [9, 15, 17]
  , [16, 18, 24]
  , [11, 17, 19]
  , [0, 18, 20]
  , [13, 19, 21]
  , [2, 20, 22]
  , [15, 21, 23]
  , [4, 22, 24]
  , [17, 23, 25]
  , [0, 6, 24]
  ]

-- | Coxeter (28 vertices, cubic symmetric, LCF [5,-5,13,-13] x 7)
fosterCoxeter :: Topology
fosterCoxeter = V.fromList
  [ [1, 5, 27]
  , [0, 2, 24]
  , [1, 3, 15]
  , [2, 4, 18]
  , [3, 5, 9]
  , [0, 4, 6]
  , [5, 7, 19]
  , [6, 8, 22]
  , [7, 9, 13]
  , [4, 8, 10]
  , [9, 11, 23]
  , [10, 12, 26]
  , [11, 13, 17]
  , [8, 12, 14]
  , [13, 15, 27]
  , [2, 14, 16]
  , [15, 17, 21]
  , [12, 16, 18]
  , [3, 17, 19]
  , [6, 18, 20]
  , [19, 21, 25]
  , [16, 20, 22]
  , [7, 21, 23]
  , [10, 22, 24]
  , [1, 23, 25]
  , [20, 24, 26]
  , [11, 25, 27]
  , [0, 14, 26]
  ]

-- | Tutte-Coxeter / Levi graph (30 vertices, cubic symmetric)
fosterTutteCoxeter :: Topology
fosterTutteCoxeter = V.fromList
  [ [1, 17, 29]
  , [0, 2, 22]
  , [1, 3, 9]
  , [2, 4, 26]
  , [3, 5, 13]
  , [4, 6, 18]
  , [5, 7, 23]
  , [6, 8, 28]
  , [7, 9, 15]
  , [2, 8, 10]
  , [9, 11, 19]
  , [10, 12, 24]
  , [11, 13, 29]
  , [4, 12, 14]
  , [13, 15, 21]
  , [8, 14, 16]
  , [15, 17, 25]
  , [0, 16, 18]
  , [5, 17, 19]
  , [10, 18, 20]
  , [19, 21, 27]
  , [14, 20, 22]
  , [1, 21, 23]
  , [6, 22, 24]
  , [11, 23, 25]
  , [16, 24, 26]
  , [3, 25, 27]
  , [20, 26, 28]
  , [7, 27, 29]
  , [0, 12, 28]
  ]

-- ---------------------------------------------------------------------------
-- Topology lookup
-- ---------------------------------------------------------------------------

buildTopology :: String -> Int -> Topology
buildTopology name n = case name of
  "disconnected"      -> disconnected n
  "ring"              -> ring n
  "star"              -> star n
  "complete"          -> complete n
  "hypercube"         -> hypercube 3  -- k=3 for 8 islands
  "barbell"           -> barbell n
  "watts-strogatz"    -> wattsStrogatz n 4 0.3 42
  "random-regular"    -> randomRegular n 3 42
  -- Foster census cubic symmetric graphs (fixed topology, ignore n)
  "k4"              -> fosterK4
  "k33"             -> fosterK33
  "cube"            -> fosterCube
  "petersen"        -> fosterPetersen
  "heawood"         -> fosterHeawood
  "mobius-kantor"   -> fosterMobiusKantor
  "pappus"          -> fosterPappus
  "dodecahedron"    -> fosterDodecahedron
  "desargues"       -> fosterDesargues
  "nauru"           -> fosterNauru
  "f26a"            -> fosterF26A
  "coxeter"         -> fosterCoxeter
  "tutte-coxeter"   -> fosterTutteCoxeter
  -- Bridge experiment: iso-spectral families (fixed topology, n=8)
  "ring-chord1"       -> ringChord1
  "ring-chord2"       -> ringChord2
  "ring-chord3"       -> ringChord3
  "star-leaf1"        -> starLeaf1
  "star-leaf2"        -> starLeaf2
  "star-leaf3"        -> starLeaf3
  -- Cycle-length experiment: two-cycle-bridge graphs (fixed topology)
  "cycle-bridge-3-3"  -> cycleBridge33
  "cycle-bridge-3-9"  -> cycleBridge39
  "cycle-bridge-9-9"  -> cycleBridge99
  -- Directed topologies (n=8, m=16, varying cycle count)
  "dag-layer"         -> dagLayer
  "dag-wide"          -> dagWide
  "lowcyc-1"          -> lowcyc1
  "bidir-ring"        -> bidirRing
  "two-cliques"       -> twoCliques
  "mesh-cyclic"       -> meshCyclic
  "dense-triangles"   -> denseTriangles
  "ring-skip2"        -> ringSkip2
  _                   -> error $ "Unknown topology: " ++ name

-- ---------------------------------------------------------------------------
-- Run and print stats
-- ---------------------------------------------------------------------------

printStats :: [Stats] -> IO ()
printStats stats = do
  putStrLn "generation,meanFitness,bestFitness,diversity"
  mapM_ (\s -> do
    putStrLn $ show (generation s)
          ++ "," ++ show (meanFitness s)
          ++ "," ++ show (bestFitness s)
          ++ "," ++ show (diversity s)
    hFlush stdout
    ) stats

-- ---------------------------------------------------------------------------
-- Parsed configuration
-- ---------------------------------------------------------------------------

data Config = Config
  { cfgDomain       :: String   -- "maze" | "onemax"
  , cfgGridSize     :: Int      -- only used for maze domain
  , cfgTopology     :: String
  , cfgNumIslands   :: Int
  , cfgPopSize      :: Int
  , cfgMigInterval  :: Int
  , cfgNumMigrants  :: Int
  , cfgTotalGens    :: Int
  , cfgSeed         :: Int
  }

-- ---------------------------------------------------------------------------
-- Main
-- ---------------------------------------------------------------------------

main :: IO ()
main = do
  args <- getArgs
  case parseArgs args of
    Just cfg -> do
      let gen  = mkStdGen (cfgSeed cfg)
          topo = buildTopology (cfgTopology cfg) (cfgNumIslands cfg)

      hPutStrLn stderr $ "Running: " ++ cfgTopology cfg
                      ++ " | domain=" ++ cfgDomain cfg
                      ++ " | islands=" ++ show (cfgNumIslands cfg)
                      ++ " | pop=" ++ show (cfgPopSize cfg)
                      ++ " | migInterval=" ++ show (cfgMigInterval cfg)
                      ++ " | migrants=" ++ show (cfgNumMigrants cfg)
                      ++ " | gens=" ++ show (cfgTotalGens cfg)
                      ++ " | seed=" ++ show (cfgSeed cfg)

      case cfgDomain cfg of
        "onemax" -> do
          let stats = runSimulation (Proxy :: Proxy OneMax)
                        (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                        (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
          printStats stats
        "maze" ->
          case cfgGridSize cfg of
            8 -> do
              let stats = runSimulation (Proxy :: Proxy Maze8)
                            (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                            (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
              printStats stats
            15 -> do
              let stats = runSimulation (Proxy :: Proxy Maze)
                            (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                            (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
              printStats stats
            gs -> do
              hPutStrLn stderr $ "Unsupported grid size: " ++ show gs ++ ". Use 8 or 15."
        "nk0" -> do
          let stats = runSimulation (Proxy :: Proxy NK0Individual)
                        (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                        (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
          printStats stats
        "nk2" -> do
          let stats = runSimulation (Proxy :: Proxy NK2Individual)
                        (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                        (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
          printStats stats
        "nk4" -> do
          let stats = runSimulation (Proxy :: Proxy NK4Individual)
                        (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                        (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
          printStats stats
        "nk6" -> do
          let stats = runSimulation (Proxy :: Proxy NK6Individual)
                        (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                        (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
          printStats stats
        other -> do
          hPutStrLn stderr $ "Unknown domain: " ++ other ++ ". Use 'maze', 'onemax', 'nk0', 'nk2', 'nk4', or 'nk6'."

      hPutStrLn stderr "Done."

    Nothing -> do
      hPutStrLn stderr "Usage: topology-sim [--domain D] [--grid N] <topology-name> <num-islands> <pop-size> <migration-interval> <num-migrants> <total-generations> <seed>"
      hPutStrLn stderr ""
      hPutStrLn stderr "Options:"
      hPutStrLn stderr "  --domain D  Domain: 'maze' (default), 'onemax', 'nk0', 'nk2', 'nk4', or 'nk6'"
      hPutStrLn stderr "  --grid N    Grid size for maze domain (default: 15, options: 8, 15)"
      hPutStrLn stderr ""
      hPutStrLn stderr "Topologies: disconnected, ring, star, complete, hypercube, barbell, watts-strogatz, random-regular"
      hPutStrLn stderr "Directed:   dag-layer, dag-wide, lowcyc-1, bidir-ring, two-cliques, mesh-cyclic, dense-triangles, ring-skip2"
      hPutStrLn stderr ""
      hPutStrLn stderr "Examples:"
      hPutStrLn stderr "  topology-sim ring 8 50 10 5 500 42                        # 15x15 maze (default)"
      hPutStrLn stderr "  topology-sim --grid 8 ring 8 50 10 5 500 42               # 8x8 maze"
      hPutStrLn stderr "  topology-sim --domain onemax ring 8 50 10 5 500 42        # OneMax"
      hPutStrLn stderr "  topology-sim --domain nk4 ring 8 50 10 5 500 42             # NK landscape (K=4)"
      hPutStrLn stderr "  topology-sim --domain onemax --grid 8 ring 8 50 10 5 500 42  # --grid ignored for onemax"

-- ---------------------------------------------------------------------------
-- Argument parsing
-- ---------------------------------------------------------------------------

-- | Parse command-line arguments, extracting optional --domain and --grid flags.
-- Flags can appear in any order before the positional arguments.
parseArgs :: [String] -> Maybe Config
parseArgs args =
  let (domain, grid, rest) = extractFlags args "maze" 15
  in case rest of
    [topoName, ni, ps, mi, nm, tg, sd] ->
      Just Config
        { cfgDomain      = domain
        , cfgGridSize    = grid
        , cfgTopology    = topoName
        , cfgNumIslands  = read ni
        , cfgPopSize     = read ps
        , cfgMigInterval = read mi
        , cfgNumMigrants = read nm
        , cfgTotalGens   = read tg
        , cfgSeed        = read sd
        }
    _ -> Nothing

-- | Extract --domain and --grid flags from args, returning defaults for unset ones.
extractFlags :: [String] -> String -> Int -> (String, Int, [String])
extractFlags ("--domain" : d : rest) _defDomain defGrid =
  extractFlags rest d defGrid
extractFlags ("--grid" : g : rest) defDomain _defGrid =
  extractFlags rest defDomain (read g)
extractFlags rest defDomain defGrid =
  (defDomain, defGrid, rest)

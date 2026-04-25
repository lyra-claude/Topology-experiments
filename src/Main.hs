module Main where

import IslandGA (Topology, Stats(..), Checkpoint(..), runSimulation, runSimulationWithGenomes)
import Maze (Maze)    -- Domain instance for 15x15 Maze
import Maze8 (Maze8)  -- Domain instance for 8x8 Maze
import OneMax (OneMax) -- Domain instance for OneMax
import NKLandscape (NK0Individual, NK2Individual, NK4Individual, NK6Individual)
import SudokuSolver (SudokuIndividual)
import Domain (Domain(..))

import qualified Data.Vector as V
import Data.Proxy (Proxy(..))
import System.Environment (getArgs)
import System.Random (mkStdGen)
import Data.Bits (xor, testBit, popCount)
import System.IO (hFlush, stdout, hPutStrLn, stderr)
import System.Directory (createDirectoryIfMissing)
import Data.List (intercalate)

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
-- Bridge experiment: iso-spectral families (constant lambda_2, varying beta_1)
-- ---------------------------------------------------------------------------

-- | Add an undirected edge between two nodes.
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
-- Interference experiment topologies (n=8, m=9, beta_1=2, kappa=2)
-- Three directed graphs with SAME node count, edge count, cycle rank,
-- and directed cycle count — but DIFFERENT cycle arrangements.
-- Tests whether non-abelian holonomy (cactus group) causes interference.
-- ---------------------------------------------------------------------------

-- | Adjacent cycles (figure-eight): two directed cycles sharing vertex 0.
-- Cycle A: 0->1->2->3->0 (length 4)
-- Cycle B: 0->4->5->6->7->0 (length 5)
-- Cycles share exactly one vertex (node 0).
interferenceAdjacent :: Topology
interferenceAdjacent = V.fromList [ [3, 7]   -- 0: receives from 3 and 7
                                  , [0]      -- 1: receives from 0
                                  , [1]      -- 2: receives from 1
                                  , [2]      -- 3: receives from 2
                                  , [0]      -- 4: receives from 0
                                  , [4]      -- 5: receives from 4
                                  , [5]      -- 6: receives from 5
                                  , [6]      -- 7: receives from 6
                                  ]

-- | Separated cycles: two directed cycles connected by a path, no shared vertices.
-- Cycle A: 0->1->2->0 (length 3, nodes {0,1,2})
-- Path: 2->3->4 (bridge)
-- Cycle B: 4->5->6->7->4 (length 4, nodes {4,5,6,7})
-- Cycles share zero vertices.
interferenceSeparated :: Topology
interferenceSeparated = V.fromList [ [2]      -- 0: receives from 2
                                   , [0]      -- 1: receives from 0
                                   , [1]      -- 2: receives from 1
                                   , [2]      -- 3: receives from 2
                                   , [3, 7]   -- 4: receives from 3 and 7
                                   , [4]      -- 5: receives from 4
                                   , [5]      -- 6: receives from 5
                                   , [6]      -- 7: receives from 6
                                   ]

-- | Nested cycles (concentric): one cycle inside another, sharing a path.
-- Outer: 0->1->2->3->4->5->6->7->0 (length 8, the full ring)
-- Shortcut: 0->4 (creates inner cycle 0->4->5->6->7->0, length 5)
-- Cycles share the path 4->5->6->7->0 (4 edges, 5 vertices).
interferenceNested :: Topology
interferenceNested = V.fromList [ [7]      -- 0: receives from 7
                                , [0]      -- 1: receives from 0
                                , [1]      -- 2: receives from 1
                                , [2]      -- 3: receives from 2
                                , [3, 0]   -- 4: receives from 3 and 0
                                , [4]      -- 5: receives from 4
                                , [5]      -- 6: receives from 5
                                , [6]      -- 7: receives from 6
                                ]

-- ---------------------------------------------------------------------------
-- Figure-eight cycle-length experiment (n=12, m=13, beta_1=2)
-- Three directed graphs varying mean cycle length while holding
-- node count, edge count, and cycle rank constant.
-- From Claudius (2026-04-25): tests cycle-length -> spectral gap.
-- ---------------------------------------------------------------------------

-- | M1: mean cycle length 4.0 (C3 + C5 at v0, pendant chain 0->7->...->11)
figEightM1 :: Topology
figEightM1 = V.fromList [ [2, 6]   -- 0: receives from 2 (C3) and 6 (C5)
                         , [0]      -- 1: receives from 0
                         , [1]      -- 2: receives from 1
                         , [0]      -- 3: receives from 0
                         , [3]      -- 4: receives from 3
                         , [4]      -- 5: receives from 4
                         , [5]      -- 6: receives from 5
                         , [0]      -- 7: receives from 0 (pendant start)
                         , [7]      -- 8: receives from 7
                         , [8]      -- 9: receives from 8
                         , [9]      -- 10: receives from 9
                         , [10]     -- 11: receives from 10
                         ]

-- | M2: mean cycle length 5.5 (C5 + C6 at v0, pendant 0->10->11)
figEightM2 :: Topology
figEightM2 = V.fromList [ [4, 9]   -- 0: receives from 4 (C5) and 9 (C6)
                         , [0]      -- 1: receives from 0
                         , [1]      -- 2: receives from 1
                         , [2]      -- 3: receives from 2
                         , [3]      -- 4: receives from 3
                         , [0]      -- 5: receives from 0
                         , [5]      -- 6: receives from 5
                         , [6]      -- 7: receives from 6
                         , [7]      -- 8: receives from 7
                         , [8]      -- 9: receives from 8
                         , [0]      -- 10: receives from 0 (pendant start)
                         , [10]     -- 11: receives from 10
                         ]

-- | M3: mean cycle length 6.5 (C7 + C6 at v0, no pendant)
figEightM3 :: Topology
figEightM3 = V.fromList [ [6, 11]  -- 0: receives from 6 (C7) and 11 (C6)
                         , [0]      -- 1: receives from 0
                         , [1]      -- 2: receives from 1
                         , [2]      -- 3: receives from 2
                         , [3]      -- 4: receives from 3
                         , [4]      -- 5: receives from 4
                         , [5]      -- 6: receives from 5
                         , [0]      -- 7: receives from 0
                         , [7]      -- 8: receives from 7
                         , [8]      -- 9: receives from 8
                         , [9]      -- 10: receives from 9
                         , [10]     -- 11: receives from 10
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
  -- Bridge experiment topologies (iso-spectral families)
  "ring-chord1"       -> ringChord1
  "ring-chord2"       -> ringChord2
  "ring-chord3"       -> ringChord3
  "star-leaf1"        -> starLeaf1
  "star-leaf2"        -> starLeaf2
  "star-leaf3"        -> starLeaf3
  -- Directed topologies (n=8, m=16, varying cycle count)
  "dag-layer"         -> dagLayer
  "dag-wide"          -> dagWide
  "lowcyc-1"          -> lowcyc1
  "bidir-ring"        -> bidirRing
  "two-cliques"       -> twoCliques
  "mesh-cyclic"       -> meshCyclic
  "dense-triangles"   -> denseTriangles
  "ring-skip2"        -> ringSkip2
  -- Interference experiment topologies (n=8, m=9, beta_1=2, kappa=2)
  "interference-adjacent"  -> interferenceAdjacent
  "interference-separated" -> interferenceSeparated
  "interference-nested"    -> interferenceNested
  -- Figure-eight cycle-length experiment (n=12, m=13, beta_1=2)
  "fig8-m1"           -> figEightM1
  "fig8-m2"           -> figEightM2
  "fig8-m3"           -> figEightM3
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
  { cfgDomain       :: String        -- "maze" | "onemax"
  , cfgGridSize     :: Int           -- only used for maze domain
  , cfgTopology     :: String
  , cfgNumIslands   :: Int
  , cfgPopSize      :: Int
  , cfgMigInterval  :: Int
  , cfgNumMigrants  :: Int
  , cfgTotalGens    :: Int
  , cfgSeed         :: Int
  , cfgDumpGenomes  :: Maybe String  -- Nothing or Just directory path
  }

-- ---------------------------------------------------------------------------
-- Main
-- ---------------------------------------------------------------------------

-- ---------------------------------------------------------------------------
-- Genome dump IO
-- ---------------------------------------------------------------------------

-- | Write genome data for one checkpoint to a CSV file.
-- Each line is one individual's genome (comma-separated locus values).
dumpGenomeFile :: FilePath -> Int -> V.Vector [Int] -> IO ()
dumpGenomeFile dir gen genomes = do
  let path = dir ++ "/gen_" ++ show gen ++ ".csv"
  writeFile path $ unlines
    [ intercalate "," (map show g)
    | g <- V.toList genomes
    , not (null g)  -- skip domains that don't support genome export
    ]

-- | Run simulation with genome dumps, writing both stats to stdout
-- and genome CSVs to the dump directory.
runAndDump :: Domain a => Proxy a -> Config -> Topology -> IO ()
runAndDump proxy cfg topo = do
  let gen = mkStdGen (cfgSeed cfg)
  case cfgDumpGenomes cfg of
    Nothing -> do
      -- Standard mode: stats only
      let stats = runSimulation proxy
                    (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                    (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
      printStats stats
    Just dumpDir -> do
      -- Genome dump mode
      let runDir = dumpDir ++ "/" ++ cfgTopology cfg ++ "_" ++ cfgDomain cfg
                   ++ "_seed" ++ show (cfgSeed cfg)
      createDirectoryIfMissing True runDir
      hPutStrLn stderr $ "Genome dump directory: " ++ runDir
      let checkpoints = runSimulationWithGenomes proxy
                          (cfgPopSize cfg) (cfgMigInterval cfg) (cfgNumMigrants cfg)
                          (cfgTotalGens cfg) topo (cfgNumIslands cfg) gen
      -- Print stats header
      putStrLn "generation,meanFitness,bestFitness,diversity"
      mapM_ (\cp -> do
        let s = cpStats cp
        putStrLn $ show (generation s)
              ++ "," ++ show (meanFitness s)
              ++ "," ++ show (bestFitness s)
              ++ "," ++ show (diversity s)
        hFlush stdout
        -- Write genome file for this generation
        dumpGenomeFile runDir (generation s) (cpGenomes cp)
        ) checkpoints

-- ---------------------------------------------------------------------------
-- Main
-- ---------------------------------------------------------------------------

main :: IO ()
main = do
  args <- getArgs
  case parseArgs args of
    Just cfg -> do
      let topo = buildTopology (cfgTopology cfg) (cfgNumIslands cfg)

      hPutStrLn stderr $ "Running: " ++ cfgTopology cfg
                      ++ " | domain=" ++ cfgDomain cfg
                      ++ " | islands=" ++ show (cfgNumIslands cfg)
                      ++ " | pop=" ++ show (cfgPopSize cfg)
                      ++ " | migInterval=" ++ show (cfgMigInterval cfg)
                      ++ " | migrants=" ++ show (cfgNumMigrants cfg)
                      ++ " | gens=" ++ show (cfgTotalGens cfg)
                      ++ " | seed=" ++ show (cfgSeed cfg)
                      ++ case cfgDumpGenomes cfg of
                           Nothing -> ""
                           Just d  -> " | dump-genomes=" ++ d

      case cfgDomain cfg of
        "onemax"  -> runAndDump (Proxy :: Proxy OneMax) cfg topo
        "nk0"     -> runAndDump (Proxy :: Proxy NK0Individual) cfg topo
        "nk2"     -> runAndDump (Proxy :: Proxy NK2Individual) cfg topo
        "nk4"     -> runAndDump (Proxy :: Proxy NK4Individual) cfg topo
        "nk6"     -> runAndDump (Proxy :: Proxy NK6Individual) cfg topo
        "maze"    ->
          case cfgGridSize cfg of
            8  -> runAndDump (Proxy :: Proxy Maze8) cfg topo
            15 -> runAndDump (Proxy :: Proxy Maze) cfg topo
            gs -> hPutStrLn stderr $ "Unsupported grid size: " ++ show gs ++ ". Use 8 or 15."
        "sudoku"  -> runAndDump (Proxy :: Proxy SudokuIndividual) cfg topo
        other     -> hPutStrLn stderr $ "Unknown domain: " ++ other
                       ++ ". Use 'maze', 'onemax', 'nk0', 'nk2', 'nk4', 'nk6', or 'sudoku'."

      hPutStrLn stderr "Done."

    Nothing -> do
      hPutStrLn stderr "Usage: topology-sim [OPTIONS] <topology> <islands> <pop> <mig-interval> <migrants> <gens> <seed>"
      hPutStrLn stderr ""
      hPutStrLn stderr "Options:"
      hPutStrLn stderr "  --domain D          Domain: 'maze' (default), 'onemax', 'nk0', 'nk2', 'nk4', 'nk6', 'sudoku'"
      hPutStrLn stderr "  --grid N            Grid size for maze domain (default: 15, options: 8, 15)"
      hPutStrLn stderr "  --dump-genomes DIR  Dump per-generation genome CSVs to DIR for PCA analysis"
      hPutStrLn stderr ""
      hPutStrLn stderr "Topologies: disconnected, ring, star, complete, hypercube, barbell, watts-strogatz, random-regular"
      hPutStrLn stderr "Directed:   dag-layer, dag-wide, lowcyc-1, bidir-ring, two-cliques, mesh-cyclic, dense-triangles, ring-skip2"
      hPutStrLn stderr ""
      hPutStrLn stderr "Examples:"
      hPutStrLn stderr "  topology-sim ring 8 50 10 5 500 42"
      hPutStrLn stderr "  topology-sim --domain nk4 ring 8 50 10 5 500 42"
      hPutStrLn stderr "  topology-sim --domain nk4 --dump-genomes results/genomes/ ring 8 50 10 5 500 42"

-- ---------------------------------------------------------------------------
-- Argument parsing
-- ---------------------------------------------------------------------------

-- | Parse command-line arguments, extracting optional flags.
parseArgs :: [String] -> Maybe Config
parseArgs args =
  let (domain, grid, dumpDir, rest) = extractFlags args "maze" 15 Nothing
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
        , cfgDumpGenomes = dumpDir
        }
    _ -> Nothing

-- | Extract flags from args, returning defaults for unset ones.
extractFlags :: [String] -> String -> Int -> Maybe String -> (String, Int, Maybe String, [String])
extractFlags ("--domain" : d : rest) _defDomain defGrid defDump =
  extractFlags rest d defGrid defDump
extractFlags ("--grid" : g : rest) defDomain _defGrid defDump =
  extractFlags rest defDomain (read g) defDump
extractFlags ("--dump-genomes" : d : rest) defDomain defGrid _defDump =
  extractFlags rest defDomain defGrid (Just d)
extractFlags rest defDomain defGrid defDump =
  (defDomain, defGrid, defDump, rest)

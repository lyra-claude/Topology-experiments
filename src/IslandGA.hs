{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}

module IslandGA
  ( Island(..)
  , Evaluated(..)
  , Topology
  , Stats(..)
  , Checkpoint(..)
  , stepIsland
  , migrate
  , computeDiversity
  , runSimulation
  , runSimulationWithGenomes
  ) where

import Domain
import Data.Proxy (Proxy(..))
import System.Random (StdGen, split, randomR)
import qualified Data.Vector as V
import Data.List (sortBy)
import Data.Ord (comparing, Down(..))

-- ---------------------------------------------------------------------------
-- Evaluated wrapper -- caches fitness alongside each individual
-- ---------------------------------------------------------------------------

-- | An individual paired with its precomputed fitness.
-- Avoids redundant fitness calls during tournament selection, migration,
-- and stats computation.
data Evaluated a = Evaluated
  { individual    :: !a
  , cachedFitness :: !Double
  } deriving (Show)

-- | Evaluate an individual and wrap it with its fitness.
evaluate :: Domain a => a -> Evaluated a
evaluate x = Evaluated x (fitness x)
{-# INLINE evaluate #-}

-- ---------------------------------------------------------------------------
-- Types
-- ---------------------------------------------------------------------------

-- | An island: population of evaluated individuals + random generator.
data Island a = Island
  { population :: !(V.Vector (Evaluated a))
  , islandGen  :: !StdGen
  }

-- | Topology as adjacency list: island index -> list of neighbor indices.
type Topology = V.Vector [Int]

-- | Statistics recorded at each checkpoint.
data Stats = Stats
  { generation   :: !Int
  , meanFitness  :: !Double
  , bestFitness  :: !Double
  , diversity    :: !Double  -- mean pairwise distance (sampled)
  } deriving (Show)

-- ---------------------------------------------------------------------------
-- Tournament selection
-- ---------------------------------------------------------------------------

-- | Tournament selection of size 3. Picks 3 random evaluated individuals,
-- returns the one with the highest cached fitness.
tournamentSelect :: V.Vector (Evaluated a) -> StdGen -> (Evaluated a, StdGen)
tournamentSelect pop gen0 =
  let n = V.length pop
      (i1, gen1) = randomR (0, n - 1) gen0
      (i2, gen2) = randomR (0, n - 1) gen1
      (i3, gen3) = randomR (0, n - 1) gen2
      e1 = pop V.! i1
      e2 = pop V.! i2
      e3 = pop V.! i3
      f1 = cachedFitness e1
      f2 = cachedFitness e2
      f3 = cachedFitness e3
      best
        | f1 >= f2 && f1 >= f3 = e1
        | f2 >= f1 && f2 >= f3 = e2
        | otherwise            = e3
  in (best, gen3)

-- ---------------------------------------------------------------------------
-- Island step (one generation)
-- ---------------------------------------------------------------------------

-- | Run one generation on a single island.
-- For each slot in the population:
--   1. Tournament-select two parents (using cached fitness)
--   2. Crossover to produce offspring
--   3. Mutate offspring
--   4. Evaluate offspring fitness once, wrap in Evaluated
stepIsland :: Domain a => Island a -> Island a
stepIsland (Island pop gen0) =
  let n = V.length pop
      (newPop, genFinal) = buildPop 0 n V.empty gen0
  in Island newPop genFinal
  where
    buildPop !i !n !acc !g
      | i >= n    = (acc, g)
      | otherwise =
          let (ep1, g1)    = tournamentSelect pop g
              (ep2, g2)    = tournamentSelect pop g1
              (child, g3)  = crossover (individual ep1) (individual ep2) g2
              (child', g4) = mutate child g3
              !evChild     = evaluate child'
          in buildPop (i + 1) n (V.snoc acc evChild) g4

-- ---------------------------------------------------------------------------
-- Migration
-- ---------------------------------------------------------------------------

-- | Migrate between islands according to topology.
--
-- For each island, select the top nMigrants individuals (by cached fitness).
-- For each neighbor of that island, replace the worst nMigrants individuals
-- in the neighbor with copies of the migrants.
migrate :: forall a. V.Vector (Island a) -> Topology -> Int -> V.Vector (Island a)
migrate islands topo nMigrants =
  let -- Precompute: for each island, the sorted population (by cached fitness, descending)
      sortedPops :: V.Vector (V.Vector (Evaluated a))
      sortedPops = V.map (\isl -> sortPopByFitness (population isl)) islands

      -- For each island, collect all incoming migrants from its neighbors
      -- Then replace the worst individuals with incoming migrants
      applyMigration :: Int -> Island a -> Island a
      applyMigration idx (Island pop g) =
        let neighbors = topo V.! idx
            -- Collect top nMigrants from each neighbor
            incoming :: [Evaluated a]
            incoming = concatMap (\nIdx ->
              let sorted = sortedPops V.! nIdx
              in V.toList (V.take nMigrants sorted)) neighbors
            nIncoming = length incoming
        in if nIncoming == 0
             then Island pop g
             else
               -- Sort current population by fitness ascending (worst first)
               let sorted = sortPopByFitnessAsc pop
                   -- Replace worst nIncoming with incoming migrants
                   nReplace = min nIncoming (V.length sorted)
                   kept = V.drop nReplace sorted
                   newPop = V.fromList (take nReplace incoming) V.++ kept
               in Island newPop g

  in V.imap applyMigration islands

-- | Sort population by cached fitness descending (best first).
sortPopByFitness :: V.Vector (Evaluated a) -> V.Vector (Evaluated a)
sortPopByFitness pop =
  let indexed = V.toList pop
      sorted = sortBy (comparing (Down . cachedFitness)) indexed
  in V.fromList sorted

-- | Sort population by cached fitness ascending (worst first).
sortPopByFitnessAsc :: V.Vector (Evaluated a) -> V.Vector (Evaluated a)
sortPopByFitnessAsc pop =
  let indexed = V.toList pop
      sorted = sortBy (comparing cachedFitness) indexed
  in V.fromList sorted

-- ---------------------------------------------------------------------------
-- Diversity
-- ---------------------------------------------------------------------------

-- | Compute diversity: mean pairwise distance on the underlying individuals.
-- If population > 20, sample 20 random pairs instead of all n*(n-1)/2.
computeDiversity :: Domain a => V.Vector (Evaluated a) -> StdGen -> (Double, StdGen)
computeDiversity pop gen0
  | n <= 1    = (0.0, gen0)
  | n <= 20   = (allPairsDiversity pop, gen0)
  | otherwise = sampleDiversity pop 20 gen0
  where
    n = V.length pop

-- | Exact mean pairwise distance for small populations.
allPairsDiversity :: Domain a => V.Vector (Evaluated a) -> Double
allPairsDiversity pop =
  let n = V.length pop
      totalPairs = n * (n - 1) `div` 2
      sumDist = sum [ distance (individual (pop V.! i)) (individual (pop V.! j))
                    | i <- [0 .. n - 2]
                    , j <- [i + 1 .. n - 1]
                    ]
  in if totalPairs == 0 then 0.0 else sumDist / fromIntegral totalPairs

-- | Sample k random pairs and compute mean distance.
sampleDiversity :: Domain a => V.Vector (Evaluated a) -> Int -> StdGen -> (Double, StdGen)
sampleDiversity pop k gen0 =
  let n = V.length pop
      (totalDist, genFinal) = sampleLoop 0 0.0 gen0
      sampleLoop !i !acc !g
        | i >= k    = (acc, g)
        | otherwise =
            let (a, g1) = randomR (0, n - 1) g
                (b, g2) = randomR (0, n - 2) g1
                -- Ensure b /= a by shifting if needed
                b' = if b >= a then b + 1 else b
                d = distance (individual (pop V.! a)) (individual (pop V.! b'))
            in sampleLoop (i + 1) (acc + d) g2
  in (totalDist / fromIntegral k, genFinal)

-- ---------------------------------------------------------------------------
-- Simulation runner
-- ---------------------------------------------------------------------------

-- | Compute stats for the entire archipelago at a given generation.
-- Uses cached fitness values -- no recomputation.
computeStats :: Domain a => Int -> V.Vector (Island a) -> StdGen -> (Stats, StdGen)
computeStats gen islands gen0 =
  let -- Gather all evaluated individuals across all islands
      allEvaluated = V.concatMap population islands
      fitnesses = V.map cachedFitness allEvaluated
      n = V.length fitnesses
      totalFit = V.foldl' (+) 0.0 fitnesses
      meanFit = if n == 0 then 0.0 else totalFit / fromIntegral n
      bestFit = if n == 0 then 0.0 else V.maximum fitnesses
      (div', gen1) = computeDiversity allEvaluated gen0
  in (Stats gen meanFit bestFit div', gen1)

-- | Run the full island-model GA simulation.
--
-- Returns a list of Stats, one per checkpoint (every migrationInterval
-- generations, plus the final generation).
runSimulation
  :: forall a. Domain a
  => Proxy a      -- ^ Proxy to fix the domain type
  -> Int          -- ^ Population size per island
  -> Int          -- ^ Migration interval (generations between migrations)
  -> Int          -- ^ Number of migrants per migration event
  -> Int          -- ^ Total generations
  -> Topology     -- ^ Island topology (adjacency list)
  -> Int          -- ^ Number of islands
  -> StdGen       -- ^ Random seed
  -> [Stats]
runSimulation _ popSize migInterval nMigrants totalGens topo numIslands gen0 =
  let -- Initialize islands
      (islands0, genInit) = initIslands numIslands popSize gen0

      -- Run generations
      loop :: Int -> V.Vector (Island a) -> StdGen -> [Stats] -> [Stats]
      loop !gen islands statsGen acc
        | gen > totalGens = reverse acc
        | otherwise =
            -- Step all islands (one generation each)
            let islands' = V.map stepIsland islands
            in if gen `mod` migInterval == 0
                 then
                   -- Checkpoint: record stats, then migrate
                   let (stats, statsGen') = computeStats gen islands' statsGen
                       islands'' = migrate islands' topo nMigrants
                   in loop (gen + 1) islands'' statsGen' (stats : acc)
                 else
                   loop (gen + 1) islands' statsGen acc

      -- Initial stats at generation 0
      (stats0, statsGen0) = computeStats 0 islands0 genInit

  in loop 1 islands0 statsGen0 [stats0]

-- | Initialize n islands, each with popSize random individuals.
-- Each individual is evaluated once and wrapped in Evaluated.
initIslands :: Domain a => Int -> Int -> StdGen -> (V.Vector (Island a), StdGen)
initIslands numIslands popSize gen0 =
  let buildIslands !i !g !acc
        | i >= numIslands = (V.fromList (reverse acc), g)
        | otherwise =
            let (pop, g') = buildPop popSize g V.empty
                (g1, g2) = split g'
            in buildIslands (i + 1) g2 (Island pop g1 : acc)
  in buildIslands 0 gen0 []

-- | Build a population of n random individuals, each evaluated and cached.
buildPop :: Domain a => Int -> StdGen -> V.Vector (Evaluated a) -> (V.Vector (Evaluated a), StdGen)
buildPop 0 g acc = (acc, g)
buildPop !n g acc =
  let (ind, g') = randomIndividual g
      !evInd = evaluate ind
  in buildPop (n - 1) g' (V.snoc acc evInd)

-- ---------------------------------------------------------------------------
-- Checkpoint with genome data (for PCA spectrum analysis)
-- ---------------------------------------------------------------------------

-- | A checkpoint that includes full genome data for PCA analysis.
data Checkpoint a = Checkpoint
  { cpStats      :: !Stats
  , cpGenomes    :: !(V.Vector [Int])  -- ^ All genomes as int lists
  } deriving (Show)

-- | Like runSimulation, but also returns genome data at each checkpoint.
-- Uses genomeToList from the Domain typeclass.
runSimulationWithGenomes
  :: forall a. Domain a
  => Proxy a -> Int -> Int -> Int -> Int -> Topology -> Int -> StdGen
  -> [Checkpoint a]
runSimulationWithGenomes _ popSize migInterval nMigrants totalGens topo numIslands gen0 =
  let (islands0, genInit) = initIslands numIslands popSize gen0

      extractGenomes :: V.Vector (Island a) -> V.Vector [Int]
      extractGenomes isls =
        V.concatMap (\isl -> V.map (genomeToList . individual) (population isl)) isls

      loop :: Int -> V.Vector (Island a) -> StdGen -> [Checkpoint a] -> [Checkpoint a]
      loop !gen islands statsGen acc
        | gen > totalGens = reverse acc
        | otherwise =
            let islands' = V.map stepIsland islands
            in if gen `mod` migInterval == 0
                 then
                   let (stats, statsGen') = computeStats gen islands' statsGen
                       genomes = extractGenomes islands'
                       islands'' = migrate islands' topo nMigrants
                   in loop (gen + 1) islands'' statsGen' (Checkpoint stats genomes : acc)
                 else
                   loop (gen + 1) islands' statsGen acc

      (stats0, statsGen0) = computeStats 0 islands0 genInit
      genomes0 = extractGenomes islands0

  in loop 1 islands0 statsGen0 [Checkpoint stats0 genomes0]

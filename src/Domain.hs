{-# LANGUAGE TypeFamilies #-}

module Domain where

import System.Random (StdGen)

-- | The Domain type class: the Kleisli interface for pluggable domains.
--
-- Each domain defines how to generate, evaluate, recombine, mutate, and
-- compare individuals. All stochastic operations thread a StdGen explicitly
-- (no MonadRandom dependency -- keeping it lightweight).
class Domain a where
  -- | Generate a random individual.
  randomIndividual :: StdGen -> (a, StdGen)

  -- | Evaluate fitness (higher is better).
  fitness :: a -> Double

  -- | Crossover two parents to produce one offspring.
  crossover :: a -> a -> StdGen -> (a, StdGen)

  -- | Mutate an individual.
  mutate :: a -> StdGen -> (a, StdGen)

  -- | Pairwise distance for diversity measurement (normalized 0-1).
  distance :: a -> a -> Double

  -- | Export genome as a list of integers (for genome dump / PCA analysis).
  -- Default: empty list (domains that don't support dumping).
  genomeToList :: a -> [Int]
  genomeToList _ = []

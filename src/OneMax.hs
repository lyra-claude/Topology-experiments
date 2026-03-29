{-# LANGUAGE BangPatterns #-}

module OneMax
  ( OneMax(..)
  , genomeLength
  ) where

import Domain

import qualified Data.Vector.Unboxed as VU
import System.Random (StdGen, randomR)

-- ---------------------------------------------------------------------------
-- Constants
-- ---------------------------------------------------------------------------

-- | Genome length for OneMax.
genomeLength :: Int
genomeLength = 100

-- ---------------------------------------------------------------------------
-- OneMax type
-- ---------------------------------------------------------------------------

-- | A OneMax genome: a binary string of length 'genomeLength'.
-- Each element is 0 or 1 (stored as Int for efficient vector ops).
newtype OneMax = OneMax { bits :: VU.Vector Int }

-- ---------------------------------------------------------------------------
-- Helpers
-- ---------------------------------------------------------------------------

-- | Build an unboxed vector of given length using a stateful generator.
-- Similar to VU.unfoldrExactN but explicit for clarity.
generateWithGen :: Int -> (StdGen -> (Int, StdGen)) -> StdGen -> (VU.Vector Int, StdGen)
generateWithGen n f gen0 =
  let go !i !g !acc
        | i >= n    = (VU.fromListN n (reverse acc), g)
        | otherwise = let (!val, !g') = f g
                      in go (i + 1) g' (val : acc)
  in go 0 gen0 []

-- ---------------------------------------------------------------------------
-- Domain instance
-- ---------------------------------------------------------------------------

instance Domain OneMax where

  -- | Generate a random binary string of length 'genomeLength'.
  randomIndividual gen0 =
    let (vec, gen') = generateWithGen genomeLength (\g -> randomR (0 :: Int, 1) g) gen0
    in (OneMax vec, gen')

  -- | Fitness = proportion of 1-bits. Optimal = 1.0 (all ones).
  fitness (OneMax v) =
    fromIntegral (VU.sum v) / fromIntegral genomeLength

  -- | Uniform crossover: each bit independently from parent 1 or parent 2
  -- with equal probability.
  crossover (OneMax p1) (OneMax p2) gen0 =
    let pickBit g =
          let (!coin, !g') = randomR (0 :: Int, 1) g
          in (coin, g')
        (coins, gen') = generateWithGen genomeLength pickBit gen0
        child = VU.zipWith3 (\c b1 b2 -> if c == 0 then b1 else b2) coins p1 p2
    in (OneMax child, gen')

  -- | Bit-flip mutation with rate 1/N per bit.
  mutate (OneMax v) gen0 =
    let mutBit g =
          let (!r, !g') = randomR (0 :: Int, genomeLength - 1) g
          in (r, g')
        (rolls, gen') = generateWithGen genomeLength mutBit gen0
        v' = VU.zipWith (\bit r -> if r == 0 then 1 - bit else bit) v rolls
    in (OneMax v', gen')

  -- | Hamming distance normalized to [0, 1].
  distance (OneMax v1) (OneMax v2) =
    let diffs = VU.sum $ VU.zipWith (\a b -> if a /= b then (1 :: Int) else 0) v1 v2
    in fromIntegral diffs / fromIntegral genomeLength

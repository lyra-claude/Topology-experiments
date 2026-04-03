{-# LANGUAGE BangPatterns #-}

module NKLandscape
  ( NK0Individual(..)
  , NK2Individual(..)
  , NK4Individual(..)
  , NK6Individual(..)
  ) where

import Domain

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector         as V
import System.Random (StdGen, mkStdGen, randomR)

-- ---------------------------------------------------------------------------
-- Constants
-- ---------------------------------------------------------------------------

-- | Genome length for NK landscape (matches OneMax).
nkGenomeLength :: Int
nkGenomeLength = 100

-- | Fixed landscape seed — separate from the GA seed.
-- Every run on the same K uses the same fitness table.
landscapeSeed :: Int
landscapeSeed = 12345

-- ---------------------------------------------------------------------------
-- Fitness table construction
-- ---------------------------------------------------------------------------

-- | A fitness table for one locus: maps a (K+1)-bit input to a fitness
-- contribution in [0,1]. Stored as a boxed vector of Doubles, length 2^(K+1).
type LocusTable = VU.Vector Double

-- | Build the complete NK fitness table: one LocusTable per locus.
-- Uses the landscapeSeed combined with the locus index for determinism.
buildFitnessTable :: Int -> V.Vector LocusTable
buildFitnessTable k =
  let tableSize = 2 ^ (k + 1)
  in V.generate nkGenomeLength $ \i ->
       let gen0 = mkStdGen (landscapeSeed * 1000003 + i * 7919)
           (tbl, _) = generateDoubles tableSize gen0
       in tbl

-- | Generate n random Doubles in [0,1] as an unboxed vector.
generateDoubles :: Int -> StdGen -> (VU.Vector Double, StdGen)
generateDoubles n gen0 =
  let go !i !g !acc
        | i >= n    = (VU.fromListN n (reverse acc), g)
        | otherwise = let (!val, !g') = randomR (0.0 :: Double, 1.0) g
                      in go (i + 1) g' (val : acc)
  in go 0 gen0 []

-- ---------------------------------------------------------------------------
-- NK fitness evaluation (shared logic)
-- ---------------------------------------------------------------------------

-- | Evaluate NK fitness for a given K value and its precomputed table.
nkFitness :: Int -> V.Vector LocusTable -> VU.Vector Int -> Double
nkFitness k table bits =
  let n = nkGenomeLength
      -- For locus i, compute the index into its table from the (K+1) bits:
      -- bits at positions i, (i+1) mod n, ..., (i+K) mod n
      locusContrib :: Int -> Double
      locusContrib i =
        let idx = foldl' (\acc j ->
                    let pos = (i + j) `mod` n
                        bit = bits `VU.unsafeIndex` pos
                    in acc * 2 + bit
                  ) 0 [0..k]
        in (table `V.unsafeIndex` i) `VU.unsafeIndex` idx
      totalFit = sum [locusContrib i | i <- [0 .. n - 1]]
  in totalFit / fromIntegral n
  where
    foldl' f z []     = z
    foldl' f z (x:xs) = let !z' = f z x in foldl' f z' xs

-- ---------------------------------------------------------------------------
-- Shared genome operations (identical to OneMax)
-- ---------------------------------------------------------------------------

-- | Build an unboxed vector of given length using a stateful generator.
generateWithGen :: Int -> (StdGen -> (Int, StdGen)) -> StdGen -> (VU.Vector Int, StdGen)
generateWithGen n f gen0 =
  let go !i !g !acc
        | i >= n    = (VU.fromListN n (reverse acc), g)
        | otherwise = let (!val, !g') = f g
                      in go (i + 1) g' (val : acc)
  in go 0 gen0 []

nkRandomIndividual :: StdGen -> (VU.Vector Int, StdGen)
nkRandomIndividual = generateWithGen nkGenomeLength (\g -> randomR (0 :: Int, 1) g)

nkCrossover :: VU.Vector Int -> VU.Vector Int -> StdGen -> (VU.Vector Int, StdGen)
nkCrossover p1 p2 gen0 =
  let (coins, gen') = generateWithGen nkGenomeLength (\g -> randomR (0 :: Int, 1) g) gen0
      child = VU.zipWith3 (\c b1 b2 -> if c == 0 then b1 else b2) coins p1 p2
  in (child, gen')

nkMutate :: VU.Vector Int -> StdGen -> (VU.Vector Int, StdGen)
nkMutate v gen0 =
  let (rolls, gen') = generateWithGen nkGenomeLength (\g -> randomR (0 :: Int, nkGenomeLength - 1) g) gen0
      v' = VU.zipWith (\bit r -> if r == 0 then 1 - bit else bit) v rolls
  in (v', gen')

nkDistance :: VU.Vector Int -> VU.Vector Int -> Double
nkDistance v1 v2 =
  let diffs = VU.sum $ VU.zipWith (\a b -> if a /= b then (1 :: Int) else 0) v1 v2
  in fromIntegral diffs / fromIntegral nkGenomeLength

-- ---------------------------------------------------------------------------
-- NK0 (K=0): each locus independent — random additive landscape (control)
-- ---------------------------------------------------------------------------

newtype NK0Individual = NK0Individual { nk0Bits :: VU.Vector Int }

{-# NOINLINE nk0Table #-}
nk0Table :: V.Vector LocusTable
nk0Table = buildFitnessTable 0

instance Domain NK0Individual where
  randomIndividual g = let (v, g') = nkRandomIndividual g in (NK0Individual v, g')
  fitness (NK0Individual v) = nkFitness 0 nk0Table v
  crossover (NK0Individual p1) (NK0Individual p2) g =
    let (c, g') = nkCrossover p1 p2 g in (NK0Individual c, g')
  mutate (NK0Individual v) g =
    let (v', g') = nkMutate v g in (NK0Individual v', g')
  distance (NK0Individual v1) (NK0Individual v2) = nkDistance v1 v2

-- ---------------------------------------------------------------------------
-- NK2 (K=2): moderate epistasis
-- ---------------------------------------------------------------------------

newtype NK2Individual = NK2Individual { nk2Bits :: VU.Vector Int }

{-# NOINLINE nk2Table #-}
nk2Table :: V.Vector LocusTable
nk2Table = buildFitnessTable 2

instance Domain NK2Individual where
  randomIndividual g = let (v, g') = nkRandomIndividual g in (NK2Individual v, g')
  fitness (NK2Individual v) = nkFitness 2 nk2Table v
  crossover (NK2Individual p1) (NK2Individual p2) g =
    let (c, g') = nkCrossover p1 p2 g in (NK2Individual c, g')
  mutate (NK2Individual v) g =
    let (v', g') = nkMutate v g in (NK2Individual v', g')
  distance (NK2Individual v1) (NK2Individual v2) = nkDistance v1 v2

-- ---------------------------------------------------------------------------
-- NK4 (K=4): high epistasis
-- ---------------------------------------------------------------------------

newtype NK4Individual = NK4Individual { nk4Bits :: VU.Vector Int }

{-# NOINLINE nk4Table #-}
nk4Table :: V.Vector LocusTable
nk4Table = buildFitnessTable 4

instance Domain NK4Individual where
  randomIndividual g = let (v, g') = nkRandomIndividual g in (NK4Individual v, g')
  fitness (NK4Individual v) = nkFitness 4 nk4Table v
  crossover (NK4Individual p1) (NK4Individual p2) g =
    let (c, g') = nkCrossover p1 p2 g in (NK4Individual c, g')
  mutate (NK4Individual v) g =
    let (v', g') = nkMutate v g in (NK4Individual v', g')
  distance (NK4Individual v1) (NK4Individual v2) = nkDistance v1 v2

-- ---------------------------------------------------------------------------
-- NK6 (K=6): very high epistasis
-- ---------------------------------------------------------------------------

newtype NK6Individual = NK6Individual { nk6Bits :: VU.Vector Int }

{-# NOINLINE nk6Table #-}
nk6Table :: V.Vector LocusTable
nk6Table = buildFitnessTable 6

instance Domain NK6Individual where
  randomIndividual g = let (v, g') = nkRandomIndividual g in (NK6Individual v, g')
  fitness (NK6Individual v) = nkFitness 6 nk6Table v
  crossover (NK6Individual p1) (NK6Individual p2) g =
    let (c, g') = nkCrossover p1 p2 g in (NK6Individual c, g')
  mutate (NK6Individual v) g =
    let (v', g') = nkMutate v g in (NK6Individual v', g')
  distance (NK6Individual v1) (NK6Individual v2) = nkDistance v1 v2

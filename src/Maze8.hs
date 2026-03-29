{-# LANGUAGE BangPatterns #-}

module Maze8
  ( Maze8(..)
  , gridSize8
  , numEdges8
  ) where

import Domain
import MazeCore

import qualified Data.Vector.Unboxed as VU
import qualified Data.IntSet as IS
import System.Random (StdGen, randomR)

-- ---------------------------------------------------------------------------
-- Constants (8x8)
-- ---------------------------------------------------------------------------

-- | Grid dimension (8x8).
gridSize8 :: Int
gridSize8 = 8

-- | Total edges in 8x8 grid: 8*7 + 7*8 = 56 + 56 = 112.
numEdges8 :: Int
numEdges8 = numEdgesFor gridSize8

-- ---------------------------------------------------------------------------
-- Maze8 type
-- ---------------------------------------------------------------------------

-- | An 8x8 maze genome: a permutation of edge indices [0..numEdges8-1].
newtype Maze8 = Maze8 { edgePermutation8 :: VU.Vector Int }

-- ---------------------------------------------------------------------------
-- Domain instance (8x8)
-- ---------------------------------------------------------------------------

instance Domain Maze8 where

  randomIndividual gen =
    let (perm, gen') = fisherYates numEdges8 gen
    in (Maze8 perm, gen')

  fitness (Maze8 perm) =
    let tree = decodeMazeFor gridSize8 perm
        pathLen = bfsShortestPathFor gridSize8 tree
    in fromIntegral pathLen / fromIntegral (numCellsFor gridSize8)

  crossover (Maze8 p1) (Maze8 p2) gen =
    let len = VU.length p1
        (i, gen1) = randomR (0, len - 2) gen
        (j, gen2) = randomR (i + 1, len - 1) gen1
        child = orderCrossover p1 p2 i j
    in (Maze8 child, gen2)

  mutate (Maze8 perm) gen =
    let len = VU.length perm
        (i, gen1) = randomR (0, len - 1) gen
        (j, gen2) = randomR (0, len - 1) gen1
        perm' = swapVec perm i j
    in (Maze8 perm', gen2)

  distance (Maze8 p1) (Maze8 p2) =
    let len = VU.length p1
        diffs = VU.sum $ VU.zipWith (\a b -> if a /= b then (1 :: Int) else 0) p1 p2
    in fromIntegral diffs / fromIntegral len

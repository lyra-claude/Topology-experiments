{-# LANGUAGE BangPatterns #-}

module MazeCore
  ( -- * Grid geometry (parameterized by grid size)
    numCellsFor
  , numHorizontalFor
  , numVerticalFor
  , numEdgesFor
  , treeEdgesFor
    -- * Edge encoding
  , edgeCellsFor
    -- * Kruskal decode
  , decodeMazeFor
    -- * Cell degrees
  , cellDegreesFor
    -- * BFS
  , bfsShortestPathFor
  , cellNeighborsFor
    -- * Shuffle
  , fisherYates
    -- * Crossover / mutation helpers
  , orderCrossover
  , swapVec
  ) where

import UnionFind (UF, makeUF, union)

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import qualified Data.IntSet as IS
import qualified Data.Sequence as Seq
import System.Random (StdGen, randomR)
import Control.Monad.ST (runST)

-- ---------------------------------------------------------------------------
-- Grid geometry (parameterized)
-- ---------------------------------------------------------------------------

numCellsFor :: Int -> Int
numCellsFor gs = gs * gs

numHorizontalFor :: Int -> Int
numHorizontalFor gs = gs * (gs - 1)

numVerticalFor :: Int -> Int
numVerticalFor gs = (gs - 1) * gs

numEdgesFor :: Int -> Int
numEdgesFor gs = numHorizontalFor gs + numVerticalFor gs

treeEdgesFor :: Int -> Int
treeEdgesFor gs = numCellsFor gs - 1

-- ---------------------------------------------------------------------------
-- Edge encoding (parameterized by grid size)
-- ---------------------------------------------------------------------------

edgeCellsFor :: Int -> Int -> (Int, Int)
edgeCellsFor !gs !e
  | e < numHorizontalFor gs =
      let r = e `div` (gs - 1)
          c = e `mod` (gs - 1)
      in (r * gs + c, r * gs + c + 1)
  | otherwise =
      let nh = numHorizontalFor gs
          offset = e - nh
          r = offset `div` gs
          c = offset `mod` gs
      in (r * gs + c, (r + 1) * gs + c)

-- ---------------------------------------------------------------------------
-- Fisher-Yates shuffle
-- ---------------------------------------------------------------------------

fisherYates :: Int -> StdGen -> (VU.Vector Int, StdGen)
fisherYates n gen0 = runST $ do
  vec <- VUM.new n
  let initLoop !i
        | i >= n    = return ()
        | otherwise = VUM.write vec i i >> initLoop (i + 1)
  initLoop 0
  gen' <- shuffleLoop vec (n - 1) gen0
  result <- VU.unsafeFreeze vec
  return (result, gen')
  where
    shuffleLoop _ 0 g = return g
    shuffleLoop vec !i g = do
      let (j, g') = randomR (0, i) g
      vi <- VUM.read vec i
      vj <- VUM.read vec j
      VUM.write vec i vj
      VUM.write vec j vi
      shuffleLoop vec (i - 1) g'

-- ---------------------------------------------------------------------------
-- Kruskal's algorithm (parameterized)
-- ---------------------------------------------------------------------------

decodeMazeFor :: Int -> VU.Vector Int -> IS.IntSet
decodeMazeFor gs perm = go 0 (makeUF (numCellsFor gs)) IS.empty 0
  where
    te = treeEdgesFor gs
    go !added !uf !edgeSet !idx
      | added >= te = edgeSet
      | idx >= VU.length perm = edgeSet
      | otherwise =
          let e = perm VU.! idx
              (cellA, cellB) = edgeCellsFor gs e
              (merged, uf') = union uf cellA cellB
          in if merged
               then go (added + 1) uf' (IS.insert e edgeSet) (idx + 1)
               else go added uf' edgeSet (idx + 1)

-- ---------------------------------------------------------------------------
-- Cell degrees (parameterized)
-- ---------------------------------------------------------------------------

-- | Compute the degree (number of passages) of each cell in a decoded maze.
-- Returns a vector indexed by cell ID.
cellDegreesFor :: Int -> IS.IntSet -> VU.Vector Int
cellDegreesFor gs treeEdges' = runST $ do
  degs <- VUM.replicate nc (0 :: Int)
  let go !idx
        | idx >= ne = return ()
        | IS.member idx treeEdges' =
            let (a, b) = edgeCellsFor gs idx
            in do VUM.modify degs (+1) a
                  VUM.modify degs (+1) b
                  go (idx + 1)
        | otherwise = go (idx + 1)
  go 0
  VU.unsafeFreeze degs
  where
    nc = numCellsFor gs
    ne = numEdgesFor gs

-- ---------------------------------------------------------------------------
-- BFS shortest path (parameterized)
-- ---------------------------------------------------------------------------

bfsShortestPathFor :: Int -> IS.IntSet -> Int
bfsShortestPathFor gs treeEdges' = bfs (Seq.singleton (0, 0)) IS.empty
  where
    nc = numCellsFor gs
    target = nc - 1

    bfs Seq.Empty _ = 0
    bfs (front Seq.:<| rest) visited
      | fst front == target = snd front
      | IS.member (fst front) visited = bfs rest visited
      | otherwise =
          let cell = fst front
              dist = snd front
              visited' = IS.insert cell visited
              neighbors = cellNeighborsFor gs cell treeEdges'
              newEntries = foldr
                (\n acc -> if IS.member n visited' then acc else acc Seq.|> (n, dist + 1))
                rest
                neighbors
          in bfs newEntries visited'

cellNeighborsFor :: Int -> Int -> IS.IntSet -> [Int]
cellNeighborsFor gs cell treeEdges' =
  let r = cell `div` gs
      c = cell `mod` gs
      gs1 = gs - 1
      nh = numHorizontalFor gs
      right = if c < gs1 && IS.member (r * gs1 + c) treeEdges'
              then [cell + 1] else []
      left  = if c > 0   && IS.member (r * gs1 + (c - 1)) treeEdges'
              then [cell - 1] else []
      down  = if r < gs1 && IS.member (nh + r * gs + c) treeEdges'
              then [cell + gs] else []
      up    = if r > 0   && IS.member (nh + (r - 1) * gs + c) treeEdges'
              then [cell - gs] else []
  in right ++ left ++ down ++ up

-- ---------------------------------------------------------------------------
-- Order Crossover (OX)
-- ---------------------------------------------------------------------------

orderCrossover :: VU.Vector Int -> VU.Vector Int -> Int -> Int -> VU.Vector Int
orderCrossover p1 p2 i j = runST $ do
  let len = VU.length p1
      segLen = j - i + 1

  let segSet = IS.fromList [p1 VU.! k | k <- [i..j]]

  let p2order = [ p2 VU.! ((j + 1 + k) `mod` len) | k <- [0 .. len - 1] ]
      fillElems = filter (\x -> not (IS.member x segSet)) p2order

  child <- VUM.new len

  let copySegment !k
        | k > j     = return ()
        | otherwise = VUM.write child k (p1 VU.! k) >> copySegment (k + 1)
  copySegment i

  let fillLoop _ [] = return ()
      fillLoop !pos (e:es)
        | pos >= len - segLen = return ()
        | otherwise = do
            let actualPos = (j + 1 + pos) `mod` len
            VUM.write child actualPos e
            fillLoop (pos + 1) es
  fillLoop 0 fillElems

  VU.unsafeFreeze child

-- ---------------------------------------------------------------------------
-- Helper
-- ---------------------------------------------------------------------------

swapVec :: VU.Vector Int -> Int -> Int -> VU.Vector Int
swapVec v i j
  | i == j    = v
  | otherwise =
      let vi = v VU.! i
          vj = v VU.! j
      in v VU.// [(i, vj), (j, vi)]

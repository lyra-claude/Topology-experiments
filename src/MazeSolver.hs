{-# LANGUAGE BangPatterns #-}

module MazeSolver
  ( MazeSolver(..)
  , mazeGridSize
  , pathLength
  , setMazeGridSize
  ) where

import Domain
import Maze (Maze(..), gridSize, numEdges, decodeMaze, bfsShortestPath)
import UnionFind ()

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import qualified Data.IntSet as IS
import System.Random (StdGen, mkStdGen, randomR)
import Control.Monad.ST (runST)

-- ---------------------------------------------------------------------------
-- Configuration
-- ---------------------------------------------------------------------------

-- | Grid size for the solver maze. Tunable difficulty parameter.
-- Larger grids = harder mazes = analogue of NK's K parameter.
-- Default: same as Maze module's gridSize (15).
--
-- The maze is deterministically generated from a fixed seed, so all
-- individuals in a run navigate the same maze.
mazeGridSize :: Int
mazeGridSize = gridSize  -- Uses the Maze module's constant (15x15)

-- | Number of steps in the path genome.
-- Longer than optimal to allow redundancy; the GA prunes inefficiency.
pathLength :: Int
pathLength = mazeGridSize * mazeGridSize * 2  -- 450 for 15x15

-- | Setter stub for future configurability.
setMazeGridSize :: Int -> IO ()
setMazeGridSize _ = return ()  -- placeholder

-- ---------------------------------------------------------------------------
-- Fixed maze (generated once from a deterministic seed)
-- ---------------------------------------------------------------------------

-- | The fixed maze that all individuals navigate. Generated via Kruskal's
-- from a deterministic seed. This is a spanning tree encoded as an IntSet
-- of edge indices.
{-# NOINLINE fixedMaze #-}
fixedMaze :: IS.IntSet
fixedMaze =
  let -- Generate a random permutation of edges using a fixed seed
      seed = mkStdGen 12345
      (perm, _) = fisherYatesMaze numEdges seed
      maze = Maze perm
  in decodeMaze maze

-- | The BFS shortest path length through the fixed maze.
{-# NOINLINE optimalPathLength #-}
optimalPathLength :: Int
optimalPathLength = bfsShortestPath fixedMaze

-- | Fisher-Yates shuffle (same as in Maze.hs but re-exported here to
-- avoid circular imports since Maze doesn't export it).
fisherYatesMaze :: Int -> StdGen -> (VU.Vector Int, StdGen)
fisherYatesMaze n gen0 = runST $ do
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
-- Direction encoding
-- ---------------------------------------------------------------------------

-- | Directions: 0=Up, 1=Down, 2=Left, 3=Right
type Direction = Int

-- | Apply a direction to a cell index, returning Nothing if the move
-- hits a wall (grid boundary) or there's no passage in the maze.
applyDirection :: Int -> Direction -> IS.IntSet -> Maybe Int
applyDirection cell dir tree =
  case dir of
    0 -> -- Up: r-1
      if r > 0 && IS.member (numH + (r - 1) * gs + c) tree
        then Just (cell - gs)
        else Nothing
    1 -> -- Down: r+1
      if r < gs1 && IS.member (numH + r * gs + c) tree
        then Just (cell + gs)
        else Nothing
    2 -> -- Left: c-1
      if c > 0 && IS.member (r * gs1 + (c - 1)) tree
        then Just (cell - 1)
        else Nothing
    3 -> -- Right: c+1
      if c < gs1 && IS.member (r * gs1 + c) tree
        then Just (cell + 1)
        else Nothing
    _ -> Nothing
  where
    gs   = mazeGridSize
    gs1  = gs - 1
    numH = gs * gs1
    r    = cell `div` gs
    c    = cell `mod` gs

-- ---------------------------------------------------------------------------
-- Genome: MazeSolver
-- ---------------------------------------------------------------------------

-- | A path genome: a sequence of directions [0..3].
newtype MazeSolver = MazeSolver { directions :: VU.Vector Int }

-- ---------------------------------------------------------------------------
-- Simulation: walk the path
-- ---------------------------------------------------------------------------

-- | Walk the path through the fixed maze, returning the final cell
-- and the number of valid (non-wall) steps taken.
walkPath :: MazeSolver -> (Int, Int)
walkPath (MazeSolver dirs) = go 0 0 0
  where
    len = VU.length dirs
    go !pos !validSteps !idx
      | idx >= len = (pos, validSteps)
      | otherwise  =
          let d = dirs VU.! idx
          in case applyDirection pos d fixedMaze of
               Just pos' -> go pos' (validSteps + 1) (idx + 1)
               Nothing   -> go pos  validSteps       (idx + 1)

-- | Manhattan distance from a cell to the target (bottom-right corner).
manhattanToExit :: Int -> Int
manhattanToExit cell =
  let gs  = mazeGridSize
      r   = cell `div` gs
      c   = cell `mod` gs
      targetR = gs - 1
      targetC = gs - 1
  in abs (r - targetR) + abs (c - targetC)

-- ---------------------------------------------------------------------------
-- Domain instance
-- ---------------------------------------------------------------------------

instance Domain MazeSolver where

  -- | Generate a random path (sequence of random directions).
  randomIndividual gen =
    let len = pathLength
        (dirs, gen') = randomDirs len gen
    in (MazeSolver dirs, gen')

  -- | Fitness: combines proximity to exit with path efficiency.
  --
  -- Components:
  --   1. Proximity (0-1): 1.0 when at exit, 0.0 at farthest corner.
  --      This is the primary signal; reaching the exit is paramount.
  --
  --   2. Efficiency bonus (0-0.5): ratio of optimal path length to
  --      actual valid steps taken. Rewards shorter paths.
  --      Only applies if the maze was solved (reached the exit).
  --
  -- Total fitness range: [0.0, 1.5]
  -- Normalized to [0.0, 1.0] by dividing by 1.5.
  fitness solver =
    let (finalCell, validSteps) = walkPath solver
        target = mazeGridSize * mazeGridSize - 1
        maxManhattan = fromIntegral (2 * (mazeGridSize - 1)) :: Double

        -- Proximity component: 1.0 at exit, 0.0 at start corner
        proximity = 1.0 - (fromIntegral (manhattanToExit finalCell) / maxManhattan)

        -- Efficiency bonus: only if solved
        efficiencyBonus
          | finalCell == target && validSteps > 0 =
              0.5 * (fromIntegral optimalPathLength / fromIntegral validSteps)
          | otherwise = 0.0

    in (proximity + efficiencyBonus) / 1.5

  -- | Two-point crossover: swap a segment of directions between parents.
  crossover (MazeSolver p1) (MazeSolver p2) gen =
    let len = VU.length p1
        (i, gen1) = randomR (0, len - 2) gen
        (j, gen2) = randomR (i + 1, len - 1) gen1
        -- Take [0..i-1] from p1, [i..j] from p2, [j+1..len-1] from p1
        child = VU.generate len (\k ->
          if k >= i && k <= j
            then p2 VU.! k
            else p1 VU.! k)
    in (MazeSolver child, gen2)

  -- | Mutate: flip a random direction to a new random direction.
  mutate (MazeSolver dirs) gen =
    let len = VU.length dirs
        (idx, gen1) = randomR (0, len - 1) gen
        (newDir, gen2) = randomR (0, 3) gen1
        dirs' = dirs VU.// [(idx, newDir)]
    in (MazeSolver dirs', gen2)

  -- | Distance: fraction of positions where directions differ.
  distance (MazeSolver d1) (MazeSolver d2) =
    let len = VU.length d1
        diffs = VU.sum $ VU.zipWith (\a b -> if a /= b then (1 :: Int) else 0) d1 d2
    in fromIntegral diffs / fromIntegral len

-- ---------------------------------------------------------------------------
-- Helper: generate random direction vector
-- ---------------------------------------------------------------------------

-- | Generate a vector of n random directions (0-3).
randomDirs :: Int -> StdGen -> (VU.Vector Int, StdGen)
randomDirs n gen0 = runST $ do
  vec <- VUM.new n
  gen' <- fillLoop vec 0 gen0
  result <- VU.unsafeFreeze vec
  return (result, gen')
  where
    fillLoop vec !i g
      | i >= n    = return g
      | otherwise = do
          let (d, g') = randomR (0 :: Int, 3) g
          VUM.write vec i d
          fillLoop vec (i + 1) g'

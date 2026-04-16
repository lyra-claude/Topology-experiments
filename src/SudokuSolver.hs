{-# LANGUAGE BangPatterns #-}

module SudokuSolver
  ( SudokuIndividual(..)
  , sudokuNumClues
  ) where

import Domain

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import System.Random (StdGen, mkStdGen, randomR)
import Control.Monad.ST (runST)

-- ---------------------------------------------------------------------------
-- Constants
-- ---------------------------------------------------------------------------

-- | Number of given clues. More clues = easier puzzle.
-- Analogous to K in NK landscapes (difficulty parameter).
-- Default: 30 for moderate difficulty.
sudokuNumClues :: Int
sudokuNumClues = 30

-- | Grid dimension.
gridDim :: Int
gridDim = 9

-- | Total cells.
totalCells :: Int
totalCells = 81

-- ---------------------------------------------------------------------------
-- Fixed puzzle (generated once from a deterministic seed)
-- ---------------------------------------------------------------------------

-- | Generate a valid complete Sudoku grid, then remove cells to create
-- a puzzle. We use a simple deterministic construction.
--
-- The fixed clues are stored as a vector of 81 Ints:
--   -1 = not a clue (evolvable position)
--   0-8 = fixed clue value (index into [1..9])
{-# NOINLINE fixedClues #-}
fixedClues :: VU.Vector Int
fixedClues =
  let -- Start with a known valid Sudoku grid (shifted-band construction)
      -- Row r, col c gets value: (3*(r `mod` 3) + r `div` 3 + c) `mod` 9
      -- This produces a valid Sudoku for any band-shift pattern.
      fullGrid = VU.generate totalCells $ \idx ->
        let r = idx `div` gridDim
            c = idx `mod` gridDim
        in (3 * (r `mod` 3) + r `div` 3 + c) `mod` gridDim

      -- Select which cells are clues using a deterministic shuffle
      -- (Fisher-Yates on cell indices, take first sudokuNumClues)
      seed = mkStdGen 54321
      (perm, _) = fisherYatesSudoku totalCells seed
      clueSet = VU.fromListN totalCells $ replicate totalCells False
      -- Mark the first sudokuNumClues positions from the permutation as clues
      cluePositions = VU.slice 0 sudokuNumClues perm

      isClue :: Int -> Bool
      isClue pos = VU.any (== pos) cluePositions

  in VU.generate totalCells $ \idx ->
       if isClue idx
         then fullGrid VU.! idx
         else (-1)

-- | The set of evolvable (non-clue) positions.
{-# NOINLINE evolvablePositions #-}
evolvablePositions :: VU.Vector Int
evolvablePositions =
  VU.fromList [ i | i <- [0 .. totalCells - 1], fixedClues VU.! i == (-1) ]

-- | Number of evolvable positions.
numEvolvable :: Int
numEvolvable = VU.length evolvablePositions

-- ---------------------------------------------------------------------------
-- Fisher-Yates shuffle
-- ---------------------------------------------------------------------------

fisherYatesSudoku :: Int -> StdGen -> (VU.Vector Int, StdGen)
fisherYatesSudoku n gen0 = runST $ do
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
-- Genome
-- ---------------------------------------------------------------------------

-- | A Sudoku individual: 81 values in [0..8] (row-major).
-- Clue positions are overwritten during fitness evaluation.
newtype SudokuIndividual = SudokuIndividual { sudokuGrid :: VU.Vector Int }

-- | Overlay fixed clues onto a genome, producing the effective grid.
overlayClues :: VU.Vector Int -> VU.Vector Int
overlayClues genome = VU.generate totalCells $ \idx ->
  let clue = fixedClues VU.! idx
  in if clue >= 0 then clue else genome VU.! idx

-- ---------------------------------------------------------------------------
-- Fitness: count constraint violations
-- ---------------------------------------------------------------------------

-- | Count duplicate violations in a group of 9 cells.
-- For each value 0-8, if it appears k times, that contributes (k-1) violations.
-- Maximum violations per group: 9 - 9 = 0 (if all unique) to 9 - 1 = 8 (all same).
groupViolations :: VU.Vector Int -> [Int] -> Int
groupViolations grid indices =
  let -- Count occurrences of each value 0-8
      counts = VU.create $ do
        c <- VUM.replicate gridDim (0 :: Int)
        mapM_ (\idx -> do
          let v = grid VU.! idx
          old <- VUM.read c v
          VUM.write c v (old + 1)
          ) indices
        return c
      -- Sum (count - 1) for counts > 1
  in VU.sum $ VU.map (\c -> max 0 (c - 1)) counts

-- | Total constraint violations across all rows, columns, and 3x3 boxes.
totalViolations :: VU.Vector Int -> Int
totalViolations grid =
  let -- Row violations
      rowV = sum [ groupViolations grid [r * gridDim + c | c <- [0..8]]
                 | r <- [0..8] ]
      -- Column violations
      colV = sum [ groupViolations grid [r * gridDim + c | r <- [0..8]]
                 | c <- [0..8] ]
      -- Box violations (3x3 boxes)
      boxV = sum [ groupViolations grid
                     [ (br * 3 + dr) * gridDim + (bc * 3 + dc)
                     | dr <- [0..2], dc <- [0..2] ]
                 | br <- [0..2], bc <- [0..2] ]
  in rowV + colV + boxV

-- ---------------------------------------------------------------------------
-- Domain instance
-- ---------------------------------------------------------------------------

instance Domain SudokuIndividual where

  -- | Generate a random grid (random values 0-8 for all 81 positions).
  randomIndividual gen =
    let (grid, gen') = generateRandom totalCells gen
    in (SudokuIndividual grid, gen')

  -- | Fitness: 1 / (1 + violations). Perfect solution = 1.0.
  fitness (SudokuIndividual genome) =
    let effective = overlayClues genome
        v = totalViolations effective
    in 1.0 / (1.0 + fromIntegral v)

  -- | Single-point crossover, respecting clue positions.
  crossover (SudokuIndividual p1) (SudokuIndividual p2) gen =
    let (point, gen') = randomR (0, totalCells - 2) gen
        child = VU.generate totalCells $ \i ->
          if i <= point then p1 VU.! i else p2 VU.! i
    in (SudokuIndividual child, gen')

  -- | Mutate: change one non-clue position to a random digit.
  mutate (SudokuIndividual genome) gen
    | numEvolvable == 0 = (SudokuIndividual genome, gen)
    | otherwise =
        let (eIdx, gen1) = randomR (0, numEvolvable - 1) gen
            pos = evolvablePositions VU.! eIdx
            (newVal, gen2) = randomR (0, gridDim - 1) gen1
            genome' = genome VU.// [(pos, newVal)]
        in (SudokuIndividual genome', gen2)

  -- | Distance: fraction of positions where values differ.
  distance (SudokuIndividual g1) (SudokuIndividual g2) =
    let diffs = VU.sum $ VU.zipWith (\a b -> if a /= b then (1 :: Int) else 0) g1 g2
    in fromIntegral diffs / fromIntegral totalCells

-- ---------------------------------------------------------------------------
-- Helper: generate random value vector
-- ---------------------------------------------------------------------------

generateRandom :: Int -> StdGen -> (VU.Vector Int, StdGen)
generateRandom n gen0 = runST $ do
  vec <- VUM.new n
  gen' <- fillLoop vec 0 gen0
  result <- VU.unsafeFreeze vec
  return (result, gen')
  where
    fillLoop vec !i g
      | i >= n    = return g
      | otherwise = do
          let (v, g') = randomR (0 :: Int, gridDim - 1) g
          VUM.write vec i v
          fillLoop vec (i + 1) g'

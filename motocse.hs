{-# OPTIONS_GHC -XBangPatterns #-}
-- Motif Occurrence Tool Comparison Suite
import System.Directory
import qualified Data.Map as M
import qualified Data.Set as S
import Data.List
import Control.Monad
import System.Environment
import System.Console.GetOpt
import qualified Data.ByteString.Char8 as BS
import qualified Data.Traversable as T
import Data.Maybe
import Debug.Trace
import Probes
import Text.ParserCombinators.Parsec

scorePredictions :: [(String, [String])] -> [(String, [String])] -> ()
scorePredictions goldStandard predictions =
    ()

flipPair (x, y) = (y, x)
pairListToAssociationSet = foldl' addPairListByFirstToSet S.empty
addPairListByFirstToSet s (f, l) = foldl' (\s' -> \v -> S.insert (f, v) s') s l

listCrossProduct x y = appendListCrossProduct [] x y
mapAndPrepend f v iv = foldl' (\x -> \y -> (f y):x ) iv v
appendListCrossProduct l [] y = l
appendListCrossProduct l (xs:x) y = (appendListCrossProduct (mapAndPrepend ((,)xs) y l) x y)

updateClassificationCounts standardSet testSet ((!tp, !fp), (!tn, !fn)) v =
    if S.member v testSet then
        if S.member v standardSet then
           (((tp + 1), fp), (tn, fn))
        else
            ((tp, (fp + 1)), (tn, fn))
    else
        if S.member v standardSet then
           ((tp, fp), (tn, (fn + 1)))
        else
            ((tp, fp), ((tn + 1), fn))

testPredictions :: [String] -> [String] -> [(String, [String])] -> [(String, [String])] -> ((Int, Int), (Int, Int))
testPredictions allTFs allProbes goldStandard compareTo =
    let
        -- In both sets of pairs, make the gene the first member, and the probe the second.
        goldStandardPairs = pairListToAssociationSet goldStandard
        compareToPairs = S.map flipPair $ pairListToAssociationSet compareTo
        allPairs = listCrossProduct allTFs allProbes
    in
      foldl' (updateClassificationCounts goldStandardPairs compareToPairs) ((0, 0), (0, 0)) allPairs

data CmdOptions = CmdOptions {
      cmdOptFSADirectory :: Maybe String,
      cmdOptInputPredictionFile :: Maybe String,
      cmdOptIncludeTFFile :: Maybe String,
      cmdOptShowHelp :: Bool
};

commandLineOptions =
    [
     Option [] ["help"] (NoArg (\opts -> opts {cmdOptShowHelp = True } ) ) "Shows this help message",
     Option [] ["prediction_file"] (ReqArg (\f -> \opts -> opts { cmdOptInputPredictionFile = Just f } ) "DIR") "The input file containing the predictions to be evaluate against the FSA gold standard",
     Option [] ["fsa_directory"] (ReqArg (\f -> \opts -> opts {cmdOptFSADirectory = Just f } ) "DIR") "The directory from which the 'gold standard' probe binding data should be read",
     Option [] ["include_tfs_file"] (ReqArg (\f opts -> opts { cmdOptIncludeTFFile = Just f }) "FILENAME") "A file containing a list of probes to use (optional)"
   ]

blankCommandLine = CmdOptions Nothing Nothing Nothing False

main = do
  cmdline <- getArgs
  case (getOpt RequireOrder commandLineOptions cmdline) of
    (_, _, errs@(_:_)) ->
        putStrLn ((showString "Errors: " . shows errs) "")
    (opts, _, _) ->
        case (foldl' (\state -> \t -> t state) blankCommandLine opts) of
          CmdOptions { cmdOptShowHelp = True } ->
              putStrLn $ usageInfo "motocse" commandLineOptions
          CmdOptions { cmdOptInputPredictionFile = Nothing } ->
              putStrLn "--prediction_file option is mandatory."
          CmdOptions { cmdOptFSADirectory = Nothing } ->
              putStrLn "--fsa_directory option is mandatory."
          CmdOptions { cmdOptFSADirectory = Just fsaDirectory, cmdOptInputPredictionFile = Just inputPredictionFile,
                       cmdOptIncludeTFFile = includeTFFFile } ->
              mainWithOptions fsaDirectory inputPredictionFile includeTFFFile

mainWithOptions fsaDirectory inputPredictionFile includeTFFile =
    do
      includeTFs <-
          T.traverse (liftM (map BS.unpack . BS.lines) . BS.readFile) includeTFFile
      files <- getDirectoryContents fsaDirectory
      let filesByGene = groupByFirst $ map (unsafeRight . parse parseOutGeneName "Filename") (filter ((/='.') . (flip (!!)) 0) files)
      let filesByGeneFiltered = maybe id (filter . flip (elem . fst)) includeTFs filesByGene
      boundByGene <-
          forM filesByGeneFiltered $ \(gene, filenames) ->
              do
                probes <- liftM concat $ mapM (loadProbesFromFile . (fsaDirectory++)) filenames
                return (gene, probes)
      predictBoundByGene <- loadGenesForProbesFromFile inputPredictionFile
      let allTFs = (nub . map fst) filesByGeneFiltered
      let allProbes = nub . concat $ map snd boundByGene
      let ((truePositives, falsePositives), (trueNegatives, falseNegatives)) =
            testPredictions allTFs allProbes boundByGene predictBoundByGene
      putStrLn $ (showString "TF count = " . shows (length allTFs) . showString " & probe count = " . shows (length allProbes)) ""
      putStrLn $ (showString "TruePositives = " . shows truePositives .
                  showString "; TrueNegatives = " . shows trueNegatives .
                  showString "; FalsePositives = " . shows falsePositives .
                  showString "; FalseNegatives = " . shows falseNegatives) ""
      putStrLn $ (showString "TPR = " . shows (100.0 * (fromIntegral truePositives) / (fromIntegral (truePositives + falseNegatives))) . showString "; FPR = " . shows (100.0 * (fromIntegral falsePositives) / (fromIntegral (falsePositives + trueNegatives)))) ""

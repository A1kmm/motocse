{-# OPTIONS_GHC -XBangPatterns #-}
-- Motif Occurrence Tool Comparison Suite
import Text.ParserCombinators.Parsec
import System.Directory
import qualified Data.Map as M
import qualified Data.Set as S
import Data.List
import Control.Monad
import Debug.Trace

fsaDirectory = "/home/andrew/Documents/fsafiles/"
inputPredictionFile = "/home/andrew/Documents/gmatim/hits.fasta"

unsafeRight (Left _) = error "Expected Left operand in unsafeEitherFirst"
unsafeRight (Right x) = x

parseOutGeneName = do
  genename <- many (noneOf "_")
  remainder <- many (anyChar)
  return (genename, genename ++ remainder)

addOrUpdateGroup v Nothing = Just [v]
addOrUpdateGroup v (Just x) = Just (v:x)
groupByFirstMap :: Ord a => [(a, b)] -> M.Map a [b]
groupByFirstMap = foldl' (\m -> \(k, v) -> M.alter (addOrUpdateGroup v) k m) M.empty

groupByFirst :: Ord a => [(a, b)] -> [(a, [b])]
groupByFirst = M.toList . groupByFirstMap

loadFromFileErrorCheck m filename =
    do
      probes <- parseFromFile m filename
      case probes of
        Left err ->
            do
              putStrLn $ (showString "Parse error: " . shows err) ""
              return []
        Right x -> return x

-- This gives use a list of strings - probes in the file.
loadProbesFromFile = loadFromFileErrorCheck (probeExtractor [])

probeExtractor prev =
    (
     do
       char '>'
       probe <- many $ noneOf " \n"
       many $ noneOf ">"
       probeExtractor (probe:prev)
    ) <|> (eof >> return prev)

scorePredictions :: [(String, [String])] -> [(String, [String])] -> ()
scorePredictions goldStandard predictions =
    ()

-- This gives us a list of pairs, the first member is the probe, and the second is the list of all genes corresponding to that probe.
loadGenesForProbesFromFile = loadFromFileErrorCheck (geneForAllProbes [])
geneForAllProbes l =
    (
     do
       char '>'
       probe <- many $ noneOf "\n"
       char '\n'
       genes <- genesForSingleProbe []
       geneForAllProbes ((probe, genes):l)
    ) <|>
    (
     do
       eof
       return l
    )
genesForSingleProbe l =
    (
     do
       geneFirst <- noneOf ">\n"
       geneRest <- many $ noneOf "\n"
       char '\n'
       genesForSingleProbe ((geneFirst:geneRest):l)
    ) <|> (return l)

flipPair (x, y) = (y, x)
pairListToAssociationSet = foldl' addPairListByFirstToSet S.empty
addPairListByFirstToSet s (f, l) = foldl' (\s' -> \v -> S.insert (f, v) s') s l

setCrossProduct x y = appendListCrossProduct [] (S.toList x) (S.toList y)
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

testPredictions :: [(String, [String])] -> [(String, [String])] -> ((Int, Int), (Int, Int))
testPredictions goldStandard compareTo =
    trace "Entering testPredictions" $ trace "Finished testPredictions" $!
    let
        -- In both sets of pairs, make the gene the first member, and the probe the second.
        goldStandardPairs = trace "Making goldStandardPairs" $ trace "Done making goldStandardPairs" $! pairListToAssociationSet goldStandard
        compareToPairs = S.map flipPair $ pairListToAssociationSet compareTo
        -- Get a set of all genes in the gold standard data.
        allGenes = S.map fst goldStandardPairs
        -- Get a set of all probes in the gold standard data.
        allProbes = S.map snd goldStandardPairs
        allPairs = setCrossProduct allGenes allProbes
    in
      foldl' (updateClassificationCounts goldStandardPairs compareToPairs) ((0, 0), (0, 0)) allPairs

main = do
  files <- getDirectoryContents fsaDirectory
  boundByGene <-
      let
          filesByGene = groupByFirst $ map (unsafeRight . parse parseOutGeneName "Filename") (filter ((/='.') . (flip (!!)) 0) files)
        in
          forM filesByGene $ \(gene, filenames) ->
              do
                probes <- liftM concat $ mapM (loadProbesFromFile . (fsaDirectory++)) filenames
                return (gene, probes)
  predictBoundByGene <-
      loadGenesForProbesFromFile inputPredictionFile
  let
      ((truePositives, falsePositives), (trueNegatives, falseNegatives)) = testPredictions boundByGene predictBoundByGene
    in
      putStrLn $ (showString "TruePositives = " . shows truePositives . showString "; TrueNegatives = " . shows trueNegatives . showString "; FalsePositives = " . shows falsePositives . showString "; FalseNegatives = " . shows falseNegatives) ""

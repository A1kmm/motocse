-- Motif Occurrence Tool Comparison Suite
import Text.ParserCombinators.Parsec
import System.Directory
import qualified Data.Map as M
import Data.List
import Control.Monad

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
loadGenesForProbesFromFile = loadFromFileErrorCheck (genesForAllProbes [])
geneForAllProbes l =
    (
     do
       char '>'
       probe <- many $ noneOf "\n"
       char '\n'
       genes <- genesForSingleProbe []
       geneForAllProbes (probe, genes):l
    ) <|>
    (
     do
       eof
       return l
    )
genesForSingleProbe l =
    (
     geneFirst <- noneOf ">\n"
     geneRest <- many $ noneOf "\n"
     char '\n'
     genesForSingleProbe ((geneFirst:geneRest):l)
    ) <|> (return l)

testPredictions goldStandard compareTo

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
  ((truePositves, falsePositives), (trueNegatives, falseNegatives)) <- testPredictions boundByGene predictBoundByGene
  putStrLn $ show boundByGene

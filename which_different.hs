import Probes
import Data.Ord
import Data.List
import System.Environment
import Control.Monad
import System.Directory
import Text.ParserCombinators.Parsec
import qualified Data.Traversable as T
import qualified Data.ByteString.Char8 as BS
import qualified Data.Set as S

canonicaliseGenesForProbe tfs probes = sortBy (comparing fst) . filter (not . null . snd) . filter (flip S.member probes . fst) .
                                  map (\(a, b) -> (a, sort (filter (flip S.member tfs) b)))

buildDiffOfGenesForProbe p1@((n1, g1):l1) p2@((n2, g2):l2) m =
  if n1 == n2
  then
      buildDiffOfGenesForProbe l1 l2 $ bilateralProbe g1 g2 (('>':n1):m)
  else if n1 < n2 then
      buildDiffOfGenesForProbe l1 p2 $ unilateralProbe '-' g1 (('>':n1):m)
  else
      buildDiffOfGenesForProbe p1 l2 $ unilateralProbe '+' g2 (('>':n2):m)

buildDiffOfGenesForProbe ((n1, g1):l1) [] m = buildDiffOfGenesForProbe l1 [] (unilateralProbe '-' g1 (('>':n1):m))
buildDiffOfGenesForProbe [] ((n1, g1):l1) m = buildDiffOfGenesForProbe [] l1 (unilateralProbe '+' g1 (('>':n1):m))
buildDiffOfGenesForProbe [] [] m = m

unilateralProbe c (a:g) m = unilateralProbe c g ((c:a):m)
unilateralProbe c [] m = m

bilateralProbe p1@(a1:g1) p2@(a2:g2) m =
    if a1 == a2
    then
        bilateralProbe g1 g2 (('=':a1):m)
    else if a1 < a2 then
         bilateralProbe g1 p2 (('-':a1):m)
    else
         bilateralProbe p1 g2 (('+':a2):m)

bilateralProbe p1@(_:_) [] m = unilateralProbe '-' p1 m
bilateralProbe [] p1@(_:_) m = unilateralProbe '+' p1 m
bilateralProbe [] [] m = m

main =
    do
      (fsaDirectory:includeTFFile:file1:file2:_) <- getArgs
      includeTFs <-
          liftM S.fromList $ liftM (map BS.unpack . BS.lines) (BS.readFile includeTFFile)
      files <- getDirectoryContents fsaDirectory
      let filesByGene = groupByFirst $ map (unsafeRight . parse parseOutGeneName "Filename") (filter ((/='.') . (flip (!!)) 0) files)
      let filesByGeneFiltered = filter (flip S.member includeTFs . fst) filesByGene
      boundByGene <-
          forM filesByGeneFiltered $ \(gene, filenames) ->
              do
                probes <- liftM concat $ mapM (loadProbesFromFile . (fsaDirectory++)) filenames
                return (gene, probes)
      let allTFs = S.fromList $ (nub . map fst) filesByGeneFiltered
      let allProbes = S.fromList $ nub . concat $ map snd boundByGene

      genes1 <- liftM (canonicaliseGenesForProbe allTFs allProbes) (loadGenesForProbesFromFile file1)
      genes2 <- liftM (canonicaliseGenesForProbe allTFs allProbes) (loadGenesForProbesFromFile file2)
      forM (reverse $ buildDiffOfGenesForProbe genes1 genes2 []) putStrLn

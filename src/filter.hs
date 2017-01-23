#!/usr/bin/runhaskell
{-# LANGUAGE Arrows #-}
-- , NoMonomorphismRestriction #-}
import           Text.XML.HXT.Arrow.ReadDocument
import           Text.XML.HXT.Core
import           System.FilePath
import           System.Posix.Files
import           System.Directory
import           System.IO
import           Control.Arrow
import           Control.Applicative
import           Control.Monad as M
import           Data.Either
import qualified Control.Exception as E
import qualified System.FilePath.Windows as W

main = do
  let f = "cfort2.vcxproj.filters"
  a <- runX $ readDocument [] f >>> g
  mapM_ m a
  where m [i,j] = do createDirectoryIfMissing True $ fst $ splitFileName n
                     -- M.when (takeExtension n== ".cpp") $ putStrLn n
                     E.handle (\e->hPutStrLn stderr $ show (e::E.SomeException))
                       $ rename i n >> createSymbolicLink n i
          where n = j </> i
       
g = proc x -> do
  y <- deep (hasName "ItemGroup") /> (isElem </ f) -< x
  i <- getAttrValue "Include" -< y
  j <- ( f /> single getText ) <<< getChildren -< y
  returnA -< convert  <$> [i, j]
  where f = hasName "Filter"

convert = joinPath . W.splitDirectories

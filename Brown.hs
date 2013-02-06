{-# LANGUAGE TemplateHaskell #-}

import Debug.Trace
import Linear
import Control.Lens
import Control.Applicative
import Control.Monad
import Data.Random
import System.Random.MWC
       
import Graphics.Gloss
import Graphics.Gloss.Interface.IO.Simulate
       
dt = 0.1
maxForce = 10000000

data SysState = SysState { _sysCoords :: [V2 Double] }
              deriving (Show)
makeLenses ''SysState

data Interaction = LennardJones { _iEpsilon, _iRadius :: Double }
                 deriving (Show)
makeLenses ''Interaction

data System = System { _sysParticleRadius  :: Double
                     , _sysBoxSize         :: V2 Double
                     , _sysGamma           :: Double
                     , _sysBeta            :: Double      -- ^ 1 / (k_B T)
                     , _sysInteraction     :: Interaction
                     }
            deriving (Show)
makeLenses ''System

defaultSystem = System { _sysParticleRadius = 20
                       , _sysBoxSize = 250
                       , _sysGamma = 5
                       , _sysBeta = 2.5
                       , _sysInteraction = LennardJones 100 20
                       }

type Energy = Double
     
interactionEnergy :: Interaction -> Double -> Energy
interactionEnergy (LennardJones eps s) r =
    -24 * eps * (2 * s^12/r^13 - s^6/r^7)

initsTails :: [a] -> [([a], a, [a])]
initsTails = go []
  where go done (x:rest) = (done, x, rest) : go (x:done) rest
        go _    []       = []

mapInitsTails :: Monad m => (([a],a,[a]) -> m a) -> [a] -> m [a]
mapInitsTails f xs = go [] xs
  where go done (x:rest) = do x' <- f (done,x,rest)
                              go (x':done) rest
        go done []       = return done

thermalize :: System -> V2 Double -> RVar (V2 Double)
thermalize sys p = do
    dx <- stdNormal
    dy <- stdNormal
    return $ p ^+^ sqrt (2 * beta * dt / gamma) *^ V2 dx dy
  where beta = sys^.sysBeta
        gamma = sys^.sysGamma

update :: System -> SysState -> RVar SysState
update sys s = do
    coords' <- mapInitsTails (\(done,s,rest)->
                        (return $ evolveParticle sys (done++rest) s)
                    >>= thermalize sys
                    >>= (return . enforcePeriodicBoundaries sys)
                    )
                    $ s^.sysCoords
    return $ set sysCoords coords' s

evolveParticle :: System -> [V2 Double] -> V2 Double -> V2 Double
evolveParticle sys others x =
    x ^+^ displ
  where f y = let dr = x ^-^ y
                  f = case interactionEnergy (sys^.sysInteraction) (norm dr) of
                         f | f > maxForce  -> maxForce
                         f | isNaN f       -> maxForce
                         f                 -> f
              in f *^ normalize dr
        displ = (dt / gamma) *^ foldl (^+^) zero (map f others)
        gamma = sys^.sysGamma

tr x = traceShow x x
   
enforcePeriodicBoundaries :: System -> V2 Double -> V2 Double
enforcePeriodicBoundaries sys x = clamp <$> sys^.sysBoxSize <*> x
  where clamp s x = max 0 $ min s x

initialConfig :: Int -> V2 Double -> SysState
initialConfig nParticles boxSize = 
    SysState { _sysCoords = [ V2 (realToFrac i * boxSize^._x / n)
                                 (realToFrac j * boxSize^._y / n)
                            | i <- [0..n], j <- [0..n]
                            ]
             }
  where n = sqrt (realToFrac nParticles)

withSystemRandomIO = withSystemRandom :: (GenIO -> IO a) -> IO a

main = do
    let c = initialConfig 64 250
    let sys = defaultSystem
    evolveInWindow sys c 10000

evolveInWindow :: System -> SysState -> Int -> IO ()
evolveInWindow sys s nFrames = withSystemRandomIO $ \mwc->
    simulateIO (InWindow "Simulation" (640,480) (0,0)) white 30
               s (return . render sys)
               (\viewport dt c->runRVar (update sys c) mwc) -- >>= \x->print x >> return x)

render :: System -> SysState -> Picture
render sys s = pictures $
    map (\p->color (dim blue)
             $ translate (realToFrac $ p^._x) (realToFrac $ p^._y)
             $ circleSolid r) (s^.sysCoords)
    ++ [line [(-rr,-rr), (sx,-rr), (sx,sy), (-rr,sy), (-rr,-rr)]]
  where rr = realToFrac $ sys^.sysParticleRadius
        V2 sx sy = fmap realToFrac (sys^.sysBoxSize) ^+^ pure rr
        r = realToFrac (sys^.sysParticleRadius) / 2

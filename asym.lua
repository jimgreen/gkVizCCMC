
--[[
  2D 10-moment simulation of asymmetric reconnection
  * Using dimensional-splitting scheme
  * Parameters are taken from
     Burch, et al. (Science, 2016)
     http://doi.org/10.1126/science.aaf2939
    to match the MMS event on Oct 16, 2015
  * A double tanh setup is used, but is slightly different from Burch et al. 2016
    and related work by Drake et al.
  * The L,M,N coordinates in the paper translate to x, -z, y in this script
--]]

----------------------
-- HELPER FUNCTIONS --
----------------------

function log(...)
   Lucee.logInfo(string.format(...))
end

function HERE()
   local info = debug.getinfo(2)
   local str = string.format("HERE: %d %s: %s",
       info.currentline, info.source, tostring(info.name))
   Lucee.logInfo(str)
end

function runUpdater(updater, tCurr, tEnd, inpList, outList, dir)
   updater:setCurrTime(tCurr)
   if inpList then
      updater:setIn(inpList)
   end
   if outList then
      updater:setOut(outList)
   end
   if dir then
      updater:setDir(dir)
   end
   return updater:advance(tEnd)
end

function runOutput(myQ, frame, tCurr, tag)
   myQ:write(string.format("%s_%d.h5", tag, frame), tCurr)
end

----------------
-- PARAMETERS --
----------------
-- characteristic parameters
              numFluids = 2
               gasGamma = 5./3.
             lightSpeed = 1.
                    mu0 = 1.
               epsilon0 = 1./mu0/(lightSpeed^2)
                  n1_n0 = 0.06
                  B1_B0 = 1.696
                     n0 = 1.
                     mi = 1.
                     qi = 1.
                  mi_me = 100.
                Ti0_Te0 = 1.374/0.105
                Ti1_Te1 = 7.73/1.288
              wpe0_wce0 = 2.5
                  beta0 = 2.958
                 Bz0_B0 = 0.099
                  pert0 = 0.1
              Bnoise_B0 = 0.0005
                  l_di0 = 1.
                  y0_Ly = -0.25
                  y1_Ly = 0.25
                     Lx = 40.96
                     Ly = 20.48
                     Nx = 1024
                     Ny = 512
-- derived parameters
                     me = mi / mi_me
                     qe = -qi
                   wpe0 = math.sqrt(n0*qe*qe/me/epsilon0)
                   wce0 = wpe0 / wpe0_wce0
                     B0 = - wce0 * me / qe
                     B1 = B0 * B1_B0
                    Bz0 = B0 * Bz0_B0
                     n1 = n0 * n1_n0
                    di0 = math.sqrt(mi/mu0/n0/qi/qi)
                    de0 = math.sqrt(me/mu0/n0/qe/qe)
                      l = l_di0 * di0
                     nc = (n0-n1) * ((B0+B1)/2.)^2 / (B1^2-B0^2)

                  pmag0 = B0*B0/2./mu0
                     p0 = pmag0*beta0
                 ptotal = pmag0 + p0
                Ttotal0 = p0/n0
                    Te0 = Ttotal0 / (1. + Ti0_Te0)
                    Ti0 = Ttotal0 - Te0
                  pmag1 = B1*B1/2./mu0
                     p1 = ptotal - pmag1
                Ttotal1 = p1/n1
                    Te1 = Ttotal1 / (1. + Ti1_Te1)
                    Ti1 = Ttotal1 - Te1

                     dB = pert0 * B0
                 Bnoise = Bnoise_B0 * B0
                   wci0 = qi*B0/mi
                    ke0 = 1./(10.*de0)
                    ki0 = 1./(10.*di0)

                     y0 = y0_Ly * Ly
                     y1 = y1_Ly * Ly
                 charge = {qe, qi}
                   mass = {me, mi}
-- domain and grid
                  lower = {-Lx/2., -Ly/2.}
                  upper = {Lx/2., Ly/2.}
                  cells = {Nx, Ny}
           periodicDirs = {0,1}
-- i/o control
                   tEnd = 40. / wci0
                 tFrame = 5. / wci0
     Lucee.IsRestarting = false
     Lucee.RestartFrame = -1
-- other computational parameters/switches
                    cfl = 0.9
                   cflm = 1.1*cfl
                limiter = "monotonized-centered"
    elcErrorSpeedFactor = 0
    mgnErrorSpeedFactor = 1
     elcErrorDampFactor = 0
     mgnErrorDampFactor = 0

-- diagnostic parameters for display

                   vAi0 = B0 / math.sqrt(mu0 * n0 * mi)
                    cs0 = math.sqrt(gasGamma*p0/n0/(mi+me))

log("%30s = %g", "mi/me", mi/me)
log("%30s = %g, %g", "Ti0, Te0", Ti0, Te0)
log("%31s = %g, %g", "Ti1, Te1", Ti1, Te1)
log("%31s = %g, %g", "p0, p1", p0, p1)
log("%30s = %g", "vAi0/c", vAi0/lightSpeed)
log("%30s = %g", "cs0/c", cs0/lightSpeed)
log("%30s = %g", "1/wci0", 1./wci0)
log("%30s = %g/wci0 = %g", "tFrame", tFrame*wci0, tFrame)
log("%30s = %g/wci0 = %g", "tEnd", tEnd*wci0, tEnd)
log("%30s = %g, %g", "Lx, Ly", Lx, Ly)
log("%30s = %g, %g", "Nx, Ny", Nx, Ny)
log("")

-----------------------
-- INITIAL CONDITION --
-----------------------

function init(x,y,z)
   local tanhys = math.tanh( (y-y0)/l ) - math.tanh( (y-y1)/l )
   local icoshs = 1. / (math.cosh((y-y0)/l))^2 - 1. / (math.cosh((y-y1)/l))^2
   local Pi = Lucee.Pi
 
   local Bxb = -B0 + 0.5 * (B1+B0) * tanhys
   -- assume perturbation vector potential has only z component, Az
   -- then dBx = dAz/dy, dBy = -dAz/dx, dBz = 0
   local Bx = Bxb + dB * math.sin(2.*Pi*x/Lx) * 2. * math.cos(2.*Pi*y/Ly) * math.sin(2.*Pi*y/Ly)
   local By = -dB * (Ly/Lx) * math.cos(2.*Pi*x/Lx) * math.sin(2.*Pi*y/Ly)^2
   local Bz = Bz0
   Bx = Bx + Bnoise*math.random()*math.random(-1,1)
   By = By + Bnoise*math.random()*math.random(-1,1)
   Bz = Bz + Bnoise*math.random()*math.random(-1,1)
 
   local T_e = Te0 + 0.5 * (Te1-Te0) * tanhys
   local T_i = Ti0 + 0.5 * (Ti1-Ti0) * tanhys

   local pmag = Bxb*Bxb/2./mu0
   local p = ptotal - pmag
   local n = p / (T_e + T_i)

   local p_e = n * T_e
   local p_i = n * T_i
   
   local rho_e = n * me
   local rho_i = n * mi

   local rhovx_e, rhovy_e = 0., 0.
   local rhovx_i, rhovy_i = 0., 0.
   local jz = -0.5 * (B1+B0) / l / mu0 * icoshs
   local jz_e = jz / (1. + T_i/T_e) -- current partition due to diamagnetic drift
   local jz_i = jz - jz_e
   local rhovz_e = jz_e * me / qe
   local rhovz_i = jz_i * mi / qi

   local Pxx_e = p_e + rhovx_e^2 / rho_e
   local Pyy_e = p_e + rhovy_e^2 / rho_e
   local Pzz_e = p_e + rhovz_e^2 / rho_e
   local Pxy_e = rhovx_e * rhovy_e / rho_e
   local Pxz_e = rhovx_e * rhovz_e / rho_e
   local Pyz_e = rhovy_e * rhovz_e / rho_e
   local Pxx_i = p_i + rhovx_i^2 / rho_i
   local Pyy_i = p_i + rhovy_i^2 / rho_i
   local Pzz_i = p_i + rhovz_i^2 / rho_i
   local Pxy_i = rhovx_i * rhovy_i / rho_i
   local Pxz_i = rhovx_i * rhovz_i / rho_i
   local Pyz_i = rhovy_i * rhovz_i / rho_i

   local Ex, Ey, Ez = 0., 0., 0.

   return rho_e, rhovx_e, rhovy_e, rhovz_e,
          Pxx_e, Pxy_e, Pxz_e, Pyy_e, Pyz_e, Pzz_e,
          rho_i, rhovx_i, rhovy_i, rhovz_i,
          Pxx_i, Pxy_i, Pxz_i, Pyy_i, Pyz_i, Pzz_i,
          Ex, Ey, Ez, Bx, By, Bz, 0., 0.
end

----------------------------
-- DECOMPOSITION AND GRID --
----------------------------
decomposition = DecompRegionCalc2D.CartGeneral {}
grid = Grid.RectCart2D {
   lower = lower,
   upper = upper,
   cells = cells,
   decomposition = decomposition,
   periodicDirs = periodicDirs,
}

----------
-- DATA --
----------
createData = function(numComponents)
   if not numComponents then
      numComponents = 28
   end
   return DataStruct.Field2D {
      onGrid = grid,
      numComponents = numComponents,
      ghost = {2, 2},
   }
end
q = createData()
qX = createData()
qNew = createData()
qDup = createData()

-- reuse array for dimensional-splitting
qY = qNew

getFields = function(myQ)
   return myQ:alias(0,10), myQ:alias(10,20), myQ:alias(20,28)
end
elc,ion,emf = getFields(q)
elcNew,ionNew,emfNew = getFields(qNew)
elcX,ionX,emfX = getFields(qX)
elcY,ionY,emfY = getFields(qY)

log("Main, large data (q, qNew, etc.) allocated...")

if applyDiff then
   qDiff = createData(1)
   if canSkipDiff then
      rhovDup = createData(3)
   end
end

if hasKField then
   keField = createData(1)
   kiField = createData(1)
   ke0 = nil
   ki0 = nil
end

log("Auxiliary data (inOutField, qDiff, etc.) allocated...")

-------------------------
-- BOUNDARY CONDITIONS --
-------------------------

function applyBc(myQ, tCurr, tEnd, dir)
   myQ:sync()
end

---------------------------------------
-- HYPERBOLIC EQUATIONS AND SOLVERS --
---------------------------------------
fluidEqn = HyperEquation.TenMoment {}
fluidEqnLax = HyperEquation.TenMoment { numericalFlux = "lax" }
emfEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor,
}

createSlvr = function(myEqn, input, output, myLimiter, updateDirections)
   local slvr = Updater.WavePropagation2D {
      onGrid = grid,
      equation = myEqn,
      -- one of no-limiter, zero, min-mod, superbee, 
      -- van-leer, monotonized-centered, beam-warming
      limiter = myLimiter,
      zeroLimiterSsBnd = zeroLimiterSsBnd,
      cfl = cfl,
      cflm = cflm,
      updateDirections = updateDirections,
   }
   slvr:setIn( {input} )
   slvr:setOut( {output} )
   return slvr
end

elcEqnSlvrX = createSlvr(fluidEqn, elc, elcX, limiter, {0})
ionEqnSlvrX = createSlvr(fluidEqn, ion, ionX, limiter, {0})
emfEqnSlvrX = createSlvr(emfEqn, emf, emfX, limiter, {0})

elcEqnSlvrY = createSlvr(fluidEqn, elcX, elcY, limiter, {1})
ionEqnSlvrY = createSlvr(fluidEqn, ionX, ionY, limiter, {1})
emfEqnSlvrY = createSlvr(emfEqn, emfX, emfY, limiter, {1})

slvrs = {
   {elcEqnSlvrX, ionEqnSlvrX, emfEqnSlvrX},
   {elcEqnSlvrY, ionEqnSlvrY, emfEqnSlvrY},
}

elcEqnLaxSlvrX = createSlvr(fluidEqnLax, elc, elcX, "zero", {0})
ionEqnLaxSlvrX = createSlvr(fluidEqnLax, ion, ionX, "zero", {0})
emfEqnLaxSlvrX = createSlvr(emfEqn, emf, emfX, "zero", {0})
              
elcEqnLaxSlvrY = createSlvr(fluidEqnLax, elcX, elcY, "zero", {1})
ionEqnLaxSlvrY = createSlvr(fluidEqnLax, ionX, ionY, "zero", {1})
emfEqnLaxSlvrY = createSlvr(emfEqn, emfX, emfY, "zero", {1})

laxSlvrs = {
   {elcEqnLaxSlvrX, ionEqnLaxSlvrX, emfEqnLaxSlvrX},
   {elcEqnLaxSlvrY, ionEqnLaxSlvrY, emfEqnLaxSlvrY},
}

qInList = {q, qX}
elcInList = {elc, elcX}
ionInList = {ion, ionX}
qOutList = {qX, qY}
elcOutList = {elcX, elcY}
ionOutList = {ionX, ionY}

----------------------------------
-- HYPERBOLIC EQUATION UPDATERS --
----------------------------------
function updateHyperEqns(tCurr, tEnd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(tEnd-tCurr)
   local useLaxFlux = false

   for dir = 0,1 do
      applyBc(qInList[dir+1], tCurr, tEnd, dir)
      for _,slvr in ipairs(slvrs[dir+1]) do
         slvr:setCurrTime(tCurr)
         local status, dtSuggested = slvr:advance(tEnd)
         myStatus = status and myStatus
         myDtSuggested = math.min(myDtSuggested, dtSuggested)
      end
      
      if (myStatus) then
         if ((fluidEqn:checkInvariantDomain(elcOutList[dir+1]) == false)
            or (fluidEqn:checkInvariantDomain(ionOutList[dir+1]) == false)) then
            log("** Negative density/pressure after updateHyperEqns along %s", tostring(dir))
            useLaxFlux = true
         end
      end

      if (not myStatus) or useLaxFlux then
         return myStatus, myDtSuggested, useLaxFlux
      end
   end

   return myStatus, myDtSuggested, useLaxFlux
end

function updateHyperEqnsLax(tCurr, tEnd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(tEnd-tCurr)

   for dir = 0,1 do
      applyBc(qInList[dir+1], tCurr, tEnd, dir)
      for _,slvr in ipairs(laxSlvrs[dir+1]) do
         slvr:setCurrTime(tCurr)
         local status, dtSuggested = slvr:advance(tEnd)
         myStatus = status and myStatus
         myDtSuggested = math.min(myDtSuggested, dtSuggested)
      end

      if not myStatus then
         return myStatus, myDtSuggested
      end
   end

   return myStatus, myDtSuggested
end

---------------------
-- SOURCE UPDATERS --
---------------------
srcSlvr = Updater.ImplicitTenMomentSrc2D {
   onGrid = grid,
   numFluids = numFluids,
   charge = charge,
   mass = mass,
   epsilon0 = epsilon0,
   linearSolver = "partialPivLu",
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor,
   elcErrorDampFactor = elcErrorDampFactor,
   mgnErrorDampFactor = mgnErrorDampFactor,
   resetNegativePressure = true,
}

-- updater to compute second-order derivative of field
diffCalc = Updater.RectSecondOrderCentralDiff2D { onGrid = grid }

function updateDiffSource(qIn, tCurr, dt)
   qIn:sync()
   diffCalc:setCurrTime(tCurr)
   diffCalc:setIn( {qIn} )
   diffCalc:setOut( {qDiff} )
   diffCalc:advance(tCurr+dt)
   if not resetDiff then
      qDiff:set(resetDiff)
   end
   qIn:accumulate(alpha*dt, qDiff)
end

-- collisioless heat-flux closure source updaters
elcCollisionlessSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = ke0,
   hasAverageWaveNumberField = hasKField,
   averageWaveNumberField = keField,
}
ionCollisionlessSrcSlvr = Updater.TenMomentLocalCollisionlessHeatFlux2D {
   onGrid = grid,
   averageWaveNumber = ki0,
   hasAverageWaveNumberField = hasKField,
   averageWaveNumberField = kiField,
}

updateKe = Updater.FieldFunction2D {
   onGrid = grid,
   inpComponents = {0},
   outComponents = {0},
   func = function(x,y,z,t,rho_e)
      local de = math.sqrt(elcMass^2/(mu0*rho_e*elcCharge^2))
      return 1/de
   end,
}
updateKe:setOut( {keField} )
updateKi = Updater.FieldFunction2D {
   onGrid = grid,
   inpComponents = {0},
   outComponents = {0},
   func = function(x,y,z,t,rho_i)
      local di = math.sqrt(ionMass^2/(mu0*rho_i*ionCharge^2))
      return 1/di
   end,
}
updateKi:setOut( {kiField} )

function updateSource(elcIn, ionIn, emfIn, tCurr, tEnd)
   local dt = tEnd - tCurr
   srcSlvr:setOut( {elcIn, ionIn, emfIn} )
   srcSlvr:setCurrTime(tCurr)
   srcSlvr:advance(tEnd)

   if hasKField then
      updateKe:setIn( {elcIn} )
      updateKe:advance(tEnd)
   end
   elcCollisionlessSrcSlvr:setOut( {elcIn} )
   elcCollisionlessSrcSlvr:setCurrTime(tCurr)
   elcCollisionlessSrcSlvr:advance(tEnd)

   if hasKField then
      updateKi:setIn( {ionIn} )
      updateKi:advance(tEnd)
   end
   ionCollisionlessSrcSlvr:setOut( {ionIn} )
   ionCollisionlessSrcSlvr:setCurrTime(tCurr)
   ionCollisionlessSrcSlvr:advance(tEnd)

   if applyDiff then
-- applying diffusion on electron momentum only
      local rhov_e = elcIn:alias(1,4)
      if (canSkipDiff) then
         rhovDup:copy(rhov_e)
      end
      local rhovx_e = elcIn:alias(1,2)
      local rhovy_e = elcIn:alias(2,3)
      local rhovz_e = elcIn:alias(3,4)
      updateDiffSource(rhovx_e, tCurr, dt)
      updateDiffSource(rhovy_e, tCurr, dt)
      updateDiffSource(rhovz_e, tCurr, dt)
      if (canSkipDiff and fluidEqn:checkInvariantDomain(elcIn) == false) then
         log(" ** Parabolic source leads to negative pressure. Will skip it.")
         rhov_e:copy(rhovDup)
      end
   end

end

----------------------------------------
-- HYPERBOLIC-EQUATION SYSTEM SOLVERS --
----------------------------------------
function updateSystem(tCurr, tEnd)
   local dthalf = 0.5*(tEnd-tCurr)

   -- check if parabolic term time-step is too large for stability
   if applyDiff then
      local cfla = alpha*dthalf/dMin^2
      if (cfla > cflm) then
         return false, dthalf*cfl/cfla, false
      end
   end

   updateSource(elc, ion, emf, tCurr, tCurr+dthalf)

   local status, dtSuggested, useLaxFlux = updateHyperEqns(tCurr, tEnd)

   updateSource(elcNew, ionNew, emfNew, tCurr, tCurr+dthalf)

   if (status and not useLaxFlux) then
      if ((fluidEqn:checkInvariantDomain(elcNew) == false)
       or (fluidEqn:checkInvariantDomain(ionNew) == false)
       or (qNew:hasNan())) then
         log("** Negative density/pressure or Nan at end of updateSystem")
         useLaxFlux = true
      end
   end

   return status, dtSuggested, useLaxFlux
end

function updateSystemLax(tCurr, tEnd)
   local dthalf = 0.5*(tEnd-tCurr)

   -- check if parabolic term time-step is too large for stability
   if applyDiff then
      local cfla = alpha*dthalf/dMin^2
      if (cfla > cflm) then
         return false, dthalf*cfl/cfla
      end
   end

   updateSource(elc, ion, emf, tCurr, tCurr+dthalf)

   local status, dtSuggested = updateHyperEqnsLax(tCurr, tEnd)

   updateSource(elcNew, ionNew, emfNew, tCurr, tCurr+dthalf)

   return status, dtSuggested
end

----------
-- MAIN --
----------
function runSimulation(tStart, tEnd, tFrame, initDt)
   local tCurr = tStart
   local frame = 1
   local tNextFrame = tCurr + tFrame
   local step = 1
   local stepInFrame = 1
   local myDt = initDt
   local dtSuggested = initDt
   local status = true
   local useLaxFlux = false
   if not outputEveryStepStart then
      outputEveryStepStart = -1
   end


   if (Lucee.IsRestarting) then
      rFileName = "q_" .. Lucee.RestartFrame .. ".h5"
      log("\nReading %s...", rFileName)
      tCurr = q:read(rFileName)
      log("Reading %s...done", rFileName)
      if not tCurr then
         tCurr = tStart + Lucee.RestartFrame * tFrame
      end
      frame = Lucee.RestartFrame + 1
      tNextFrame = tCurr + tFrame
      log('\nRestarting from frame %d tCurr = %g\n', Lucee.RestartFrame, tCurr)
   else
      log("Initializing...")
      q:set(init)
      log("Initializing...done")
   end
   q:sync()
   if (not Lucee.IsRestarting) then
      log(">>> Writing output 0 at t = %g...", tCurr)
      runOutput(q, 0, tStart, "q")
   end
   applyBc(q, 0, 0, 0)
   applyBc(q, 0, 0, 1)

   log("Time-integration loop starting...")

   while true do
      if alwaysUseLaxFlux then
         useLaxFlux = true
      end

      qDup:copy(q)

      if (tCurr + myDt > tEnd) then
         myDt = tEnd - tCurr
      end

      if useLaxFlux then
         log("Taking step %5d (%4d in frame %3d); t = %10g, dt = %10g; using Lax fluxes",
               step, stepInFrame, frame, tCurr, myDt)
         status, dtSuggested = updateSystemLax(tCurr, tCurr+myDt)
         if status then
            useLaxFlux = false
         end
      else
         log("Taking step %5d (%4d in frame %3d); t = %10g, dt = %10g",
               step, stepInFrame, frame, tCurr, myDt)
         status, dtSuggested, useLaxFlux = updateSystem(tCurr, tCurr+myDt)
      end

      if not status then
         log (" ** dt %g too large! Will retake step with dt %g; useLaxFlux = %s",
                            myDt, dtSuggested, tostring(useLaxFlux))
         myDt = dtSuggested
         q:copy(qDup)
      elseif useLaxFlux then
         log(" ** Negative pressure or density! Will retake step with Lax fluxes")
         q:copy(qDup)
      else

         local badElc = fluidEqn:checkInvariantDomain(elcNew) == false
         local badIon = fluidEqn:checkInvariantDomain(ionNew) == false
         local badQ = qNew:hasNan()
  
         if (badElc or badIon or badQ) then
            log(" ** Error occurred!")
            if badElc then log(" ** Negative electron pressure/density!") end
            if badIon then log(" ** Negative ion pressure/density!") end
            if badQ then log(" ** Nans occurred!") end
            log(" ** Stopping simulation, dumping this and last states!")
            runOutput(qDup, -1, tCurr, "dead")
            runOutput(qNew, 0, tCurr+myDt, "dead")
            break
         end

         q:copy(qNew)
         tCurr = tCurr + myDt

         if (outputEveryStepStart > 0 and step > outputEveryStepStart) then
            runOutput(qNew, step, tCurr, "step")
         end

         if (tCurr > tNextFrame or tCurr >= tEnd) then
            log(">>> Writing output %d at t = %g...", frame, tCurr)
            runOutput(qNew, frame, tCurr, "q")
            frame = frame + 1
            tNextFrame = tNextFrame + tFrame
            stepInFrame = 0
         end

         myDt = dtSuggested
         step = step + 1
         stepInFrame = stepInFrame + 1

         if (tCurr >= tEnd) then
            break
         end
      end
   end
end

------------------------
-- RUN THE SIMULATION --
------------------------
tStart = 0
initDt = 100
runSimulation(tStart, tEnd, tFrame, initDt)

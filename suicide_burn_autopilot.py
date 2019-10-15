# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Suicide Burn Autopilot
# 
# Copyright (c) 2019 Marc Seifert
# This Python script is published under the MIT license.
# The license text can be found at the end of this file.
#
# "Suicide Burn Autopilot" is a Python script using the 
# Kerbal Space Program mods KRPC and MechJeb to land a
# spacecraft by performing a suicide burn maneuver.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import krpc
import time
import numpy as np
from __future__ import division
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize


# ---
# constants and settings
# ---

g = 9.81335163116455    # standard gravity in KSP
MAX_STEP_IVP = 0.25     # max step size of solve_ivp in seconds
FATOL = 0.5             # errors acceptable for the convergence of minimize
XATOL = 0.05
RPC_PORT = 50002        # ports used for the KRPC connection
STREAM_PORT = 50003     # default: RPC_PORT = 50000, STREAM_PORT = 50001
CONSOLE_OUTPUT = True   # enables text output
FINAL_ALTITUDE = 10     # altitude offset from surface in m


# ---
# functions
# ---

def delta_h(t, orbit, body, refframe):
    '''
    Calculates the altitude difference
    between orbit and surface at time t.
    '''
    pos = orbit.position_at(t, refframe)
    lat = body.latitude_at_position(pos, refframe)
    lon = body.longitude_at_position(pos, refframe)
    surface_h = body.surface_height(lat, lon)
    dh = body.altitude_at_position(pos, refframe) - surface_h
    return dh


def EOM_3d(t, w, c, t_burn, refframes, sc):
    '''
    Defines the 3d equations of motion of a rocket in a gravitational field
    burning retrograde in the rotating frame of the body.
    '''
    #x, y, z, vx, vy, vz = w    # initial coordinates and velocity (rotated frame)
    #mu, F, m0, dm = c          # graviational parameter and rocket-defining parameters
    #body_fixed_frame, refframe_rotating = refframes   # body reference frames

    # convert velocity to rotating frame
    vr = sc.transform_velocity(w[0:3], w[3:6], refframes[0], refframes[1])
    
    # normalize velocity vector
    vr = np.array(vr)/np.linalg.norm(vr)
    
    # transform velocity direction to fixed frame -> retrograde direction
    rg = sc.transform_direction(vr, refframes[1], refframes[0])
    
    # distance between vessel and center of body
    d = np.linalg.norm(w[0:3])
    
    # first order equations of motion
    if t < t_burn:      # coast until t_burn
        f = [w[3],
             w[4],
             w[5],
             -c[0]*w[0]/d**3,
             -c[0]*w[1]/d**3,
             -c[0]*w[2]/d**3]
    else:               # actual hoverslam maneuver
        m = c[2]-c[3]*(t-t_burn)
        f = [w[3],
             w[4],
             w[5],
             -c[1]*rg[0]/m - c[0]*w[0]/d**3,
             -c[1]*rg[1]/m - c[0]*w[1]/d**3,
             -c[1]*rg[2]/m - c[0]*w[2]/d**3]
    return f


def cost_function(t_burn, w, constants, t_span, refframes, sc, body):
    '''
    Calculates the final hight and velocity after hoverslam maneuver.
    '''
    mu, F, m0, dm = constants
    refframe_fixed, refframe_rotating = refframes   # body reference frames
    
    # The solver evaluates EOM() within the time interval t_span. If the velocity
    # drops to zero (calculated by finish_burn) or the vessel hits the surface
    # (guess which function...), the solver stops the evaluation.
    s = solve_ivp(lambda t, w: EOM_3d(t, w, constants, t_burn, refframes, sc),
                  t_span, w, max_step=MAX_STEP_IVP, events=(finish_burn, hit_surface),
                  dense_output=True, method='LSODA')

    # only one of both events can happen
    for i in xrange(len(s.t_events)):
        if len(s.t_events[i]) > 0:
            # t_event is the time of the end of the maneuver
            t_event = s.t_events[i][0]
            break
    
    # position and velocity at impact or zero velocity
    f = s.sol(t_event)
    
    # convert velocity to rotating frame
    pos = f[0:3]
    vr = sc.transform_velocity(pos, f[3:6], refframe_fixed, refframe_rotating)
    
    # current altitude above terrain at t_event
    lat = body.latitude_at_position(pos, refframe_fixed)        # latitude
    lon = body.longitude_at_position(pos, refframe_fixed)       # longitude
    body_rotation = t_event*body.rotational_speed*180./np.pi
    lon -= body_rotation    # compensate body rotation
    if lon < -180.:
        lon += 360.
    current_alt = body.altitude_at_position(pos, refframe_fixed) - body.surface_height(lat, lon)
    
    # figure of merit: final height + final velocity relative to the surface
    FOM = 0.2*current_alt + np.linalg.norm(vr)
    
    # worse t_burn should give worse FOM (makes the minimization problem convex)
    if t_event < t_burn:
        FOM += 30*(t_burn-t_event)
    
    if CONSOLE_OUTPUT:
        print '%4.2f \t %.1f \t %.1f \t\t %.1f' % (t_burn[0], FOM, np.linalg.norm(vr), current_alt)
    return FOM


# -------
# initialize connection to KSP
conn = krpc.connect(name='Suicide Burn Autopilot', rpc_port=RPC_PORT, stream_port=STREAM_PORT)

# ---
# set up UI
# ---
canvas = conn.ui.stock_canvas
screen_size = canvas.rect_transform.size
panel = canvas.add_panel()
panel_size = (200, 120)
panel.rect_transform.size = panel_size
panel.rect_transform.position = (screen_size[0]//2-panel_size[0]//2-10, 0)
title_panel = panel.add_panel()
title_panel.rect_transform.size = (panel_size[0], 24)
title_panel.rect_transform.position = (0, panel_size[1]//2-12)
title = title_panel.add_text('Suicide Burn Autopilot UI')
title.alignment = conn.ui.TextAnchor.middle_center
title.color = (1.0, 1.0, 1.0)
title.size = 13

# add a button to start the hoverslam procedure
button = panel.add_button('Start autopilot')
button.rect_transform.size = (panel_size[0]-20, 30)
button.rect_transform.position = (0, 12)

# add stream to monitor button state
button_clicked = conn.add_stream(getattr, button, 'clicked')

# add text output
info_text = panel.add_text('')
info_text.rect_transform.position = (0, -40)
info_text.rect_transform.size = (panel_size[0]-20, 60)
info_text.alignment = conn.ui.TextAnchor.upper_left
info_text.color = (1.0, 1.0, 1.0)
info_text.size = 13


# ---
# main loop
# ---
while True:
    if button_clicked():
        button.clicked = False  # reset button state
        time.sleep(0.1)
        
        # collect information about vessel and orbit
        # cancel action if requirements for hoverslam are not fulfilled
        sc = conn.space_center
        vessel = sc.active_vessel
        orbit = vessel.orbit
        body = orbit.body
        if vessel.situation != sc.VesselSituation.sub_orbital:
            info_text.content = 'Spacecraft status not sub-orbital!'
            info_text.color = (1.0, 0.501, 0.078)
        elif body.has_atmosphere:
            info_text.content = 'Planet/moon has an atmosphere!'
            info_text.color = (1.0, 0.501, 0.078)
        elif vessel.available_thrust < 1e-6:
            info_text.content = 'No thrust available!'
            info_text.color = (1.0, 0.501, 0.078)
        else:
            # conditions checked
            # collect further information about vessel and body
            body_fixed_frame = body.non_rotating_reference_frame
            body_rotating_frame = body.reference_frame
            mu = body.gravitational_parameter
            m0 = vessel.mass
            F = vessel.max_vacuum_thrust
            dm = F/(vessel.vacuum_specific_impulse*g)    # fuel consumption (kg/s)
            
            # data streams (positions, velocities and time)
            pos = conn.add_stream(vessel.position, body_fixed_frame)
            pos_rot = conn.add_stream(vessel.position, body_rotating_frame)
            vel_fixed = conn.add_stream(vessel.velocity, body_fixed_frame)
            vel_rotating = conn.add_stream(vessel.velocity, body_rotating_frame)
            ut = conn.add_stream(getattr, sc, 'ut')    # time
            
            # current state used for calculations
            t0 = ut()
            pos0 = pos()
            vel0 = vel_fixed()
            
            # starting calculations
            info_text.content = 'Calculating hoverslam trajectory...'
            info_text.color = (1.0, 1.0, 1.0)
            
            # find the time of impact on the surface by finding the root of
            # the function delta_h using scipy.optimize.brentq
            # interval [t0, t1] in which delta_h changes its sign
            t1 = t0 + orbit.time_to_periapsis
            t_impact = brentq(delta_h, args=(orbit, body, body_fixed_frame),
                              a=t0, b=t1)
            
            # estimate time until burn
            # the factor 0.8 is empirical
            t_burn_guess = t_impact - 0.8*np.linalg.norm(vel_rotating())/(F/m0)
            
            # ---
            # NOTE: Times until now are expressed in KSP universal time (ut). For
            # the calculations it is easier to proceed in relative time starting
            # at t0 (see above). To shift to the relative time frame, t0 is 
            # subtracted from ut.
            # ---
            
            # define event functions for the optimization algorithm
            # done here to enable access to objects like sc, body and so on
            def finish_burn(t, y):
                '''
                Calculates the current speed in the rotating frame. Terminates
                the solver when reaching 1 m/s. 0 m/s is not possible as the
                function has to change sign.
                '''
                vr = sc.transform_velocity(y[0:3], y[3:6],
                                           body_fixed_frame, body_rotating_frame)
                return np.linalg.norm(vr) - 1.0
            finish_burn.terminal = True
            finish_burn.direction = -1
            
            def hit_surface(t, y):
                '''
                Calculates the current altitude. Terminates the solver when
                reaching FINAL_ALTITUDE m.
                '''
                lat = body.latitude_at_position(y[0:3], body_fixed_frame)
                lon = body.longitude_at_position(y[0:3], body_fixed_frame)
                body_rotation = t*body.rotational_speed*180./np.pi
                lon -= body_rotation
                if lon < -180.:
                    lon += 360.
                surface_h = body.surface_height(lat, lon)
                return body.altitude_at_position(y[0:3], body_fixed_frame) - surface_h - FINAL_ALTITUDE
            hit_surface.terminal = True
            hit_surface.direction = -1
            
            
            # ---
            # the actual minimization process
            # ---
            
            if CONSOLE_OUTPUT:
                print '\nminimizing cost function\n'
                print 't_burn \t FOM \t final speed \t final altitude'
            
            # w0: initial state of the vessel
            w0 = pos0 + vel0
            
            # c: constants required for optimization
            c = (mu, F, m0, dm)
            
            # t_span: time interval in which the EOM is evaluated
            t_span = (0.0, orbit.time_to_periapsis)
            
            refframes = (body_fixed_frame, body_rotating_frame)
            tt = time.time()    # measure the time of the optimization
            
            # minimize cost_function (and therefore the figure of merit "speed + altitude")
            # by varying t_burn
            res = minimize(cost_function, t_burn_guess-t0,
                           args=(w0, c, t_span, refframes, sc, body),
                           method='Nelder-Mead',
                           options={'fatol': FATOL, 'xatol': XATOL})

            # result of the minimization
            t_burn = res.x + t0
            if type(t_burn) == np.ndarray:
                t_burn = float(t_burn)
            
            # solve EOM again using the optimized t_burn
            s = solve_ivp(lambda t, w: EOM_3d(t, w, c, t_burn-t0, refframes, sc),
                          t_span, w0, max_step=MAX_STEP_IVP, events=(finish_burn, hit_surface),
                          dense_output=True, method='LSODA')
            for i in xrange(len(s.t_events)):
                if len(s.t_events[i]) > 0:
                    t_touchdown = s.t_events[i][0]
                    break
            
            # evaluate solution at end of burn (t_touchdown)
            final = s.sol(t_touchdown)  # pos and vel of vessel at end of maneuver
            posf = final[0:3]
            lat = body.latitude_at_position(posf, body_fixed_frame)
            lon = body.longitude_at_position(posf, body_fixed_frame)
            body_rotation = t_touchdown*body.rotational_speed*180./np.pi
            lon -= body_rotation    # compensate rotation of the body
            if lon < -180.0:
                lon += 360.0
            surface_h = body.surface_height(lat, lon)
            hf = body.altitude_at_position(posf, body_fixed_frame) - surface_h
            
            # final velocity
            vrf = sc.transform_velocity(posf, final[3:6],
                                        body_fixed_frame, body_rotating_frame)
            
            # required delta_v, calculated using the rocket equation
            delta_v = vessel.vacuum_specific_impulse*g*np.log(m0/(m0-dm*(t_touchdown-t_burn+t0)))
            
            if CONSOLE_OUTPUT:
                print '\nOptimization finished'
                print 'Elapsed time: %3.1f s' % (time.time()-tt)
                print '\n\nResults of the optimization:'
                print '\nt_burn = %4.2f s' % (t_burn-t0)
                print 't_touchdown = %4.2f s' % t_touchdown
                print '\nfinal altitude: %.1f m' % hf
                print 'final velocities: vx=%2.1f m/s, vy=%2.1f m/s, vz=%2.1f m/s' % (vrf[0], vrf[1], vrf[2])
                print 'final coordinates: (%3.4f, %3.4f)\n' % (lat, lon)
                print 'required delta_v: %3.1f m/s\n' % delta_v
            
            info_text.content = 'Calculation finished'
            info_text.color = (1.0, 1.0, 1.0)
            
            # if game is paused, wait
            while conn.krpc.paused:
                time.sleep(0.2)
            
            # initialize krpc and mechjeb autopilots
            ap = vessel.auto_pilot
            ap.sas = False
            ap.reference_frame = body_rotating_frame
            control = vessel.control
            mj = conn.mech_jeb
            
            info_text.content = 'Starting autopilot'
            info_text.color = (1.0, 1.0, 1.0)
            
            # approximate direction of velocity at t_burn
            pos1 = np.array(orbit.position_at(t_burn, body_rotating_frame))
            pos2 = np.array(orbit.position_at(t_burn+1.0, body_rotating_frame))
            dpos = pos1 - pos2
            v_burn = dpos/np.linalg.norm(dpos)
            
            # rotate vessel to retrograde direction at t_burn
            ap.target_direction = v_burn
            ap.engage()
            ap.wait()
            ap.disengage()
            ap.sas = True
            ap.sas_mode = ap.sas_mode.stability_assist
            
            # warp to t_burn-10s
            info_text.content = 'Warping to burn'
            info_text.color = (1.0, 1.0, 1.0)
            sc.warp_to(t_burn - 10.0)
            
            # point vessel retrograde during burn
            ap.sas_mode = ap.sas_mode.retrograde
            
            # wait until t_burn
            while ut() < t_burn:
                time.sleep(0.01)
            
            # light the candle!
            if CONSOLE_OUTPUT:
                print '\nSuicide burn started'
            control.throttle = 1.0
            control.gear = True
            info_text.content = 'Suicide burn started'
            info_text.color = (1.0, 1.0, 1.0)
            
            # wait until the vessel is slowed down below 3 m/s
            while np.linalg.norm(vel_rotating()) > 3.0:
                time.sleep(0.01)
            control.throttle = 0.0
            if CONSOLE_OUTPUT:
                print '\nSuicide burn finished'
            
            # switch to mechjeb landing autopilot
            if CONSOLE_OUTPUT:
                print '\nSwitching to MechJeb landing autopilot'
            lap = mj.landing_autopilot
            lap.enabled = True
            lap.touchdown_speed = 1.0
            lap.land_untargeted()
            
            # cut throttle and enable SAS when landed
            while vessel.situation != sc.VesselSituation.landed:
                time.sleep(0.01)
            control.throttle = 0.0
            lap.enabled = False
            ap.sas = True
            ap.sas_mode = ap.sas_mode.stability_assist
            
            info_text.content = 'Somehow we landed'
            info_text.color = (0.0, 1.0, 0.0)
            
            if CONSOLE_OUTPUT:
                print 'landed'
            
            # the actual landing coordinates
            lat = body.latitude_at_position(pos(), body_fixed_frame)
            lon = body.longitude_at_position(pos(), body_fixed_frame)
            if CONSOLE_OUTPUT:
                print '\nlanding coordinates: (%3.4f, %3.4f)\n' % (lat, lon)
            
            # wait until situation changes
            while vessel.situation == sc.VesselSituation.landed:
                time.sleep(1)
            info_text.content = ''


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# MIT License
#
# Copyright (c) 2019 Marc Seifert
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial
# portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #